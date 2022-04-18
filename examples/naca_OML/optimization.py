#!/usr/bin/env python
"""
This file is part of the package FUNtoFEM for coupled aeroelastic simulation
and design optimization.

Copyright (C) 2015 Georgia Tech Research Corporation.
Additional copyright (C) 2015 Kevin Jacobson, Jan Kiviaho and Graeme Kennedy.
All rights reserved.

FUNtoFEM is licensed under the Apache License, Version 2.0 (the "License");
you may not use this software except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

#FUNtoFEM software is used here for full aerothermoelastic optimization of a NACA wing
#   with parametric geometries built by ESP/CAPS
#Authors: Sean Engelstad, Sejal Sahu, Graeme Kennedy
#Work for Aviation Paper 2022 on Full aerothermoelastic optimization with Parametric geometries of ESP/CAPS

from __future__ import print_function

#import funtofem classes and functions
from pyfuntofem.model  import *
from pyfuntofem.driver import *
from pyfuntofem.fun3d_interface import *
from pyfuntofem.tacs_interface import TacsSteadyInterface

#import from tacs
from tacs.pytacs import pyTACS
from tacs import functions
from tacs import TACS, functions, constitutive, elements, pyTACS, problems

#import normal classes
import os, shutil
import numpy as np

#import other modules
from pyoptsparse import SLSQP, Optimization
import pyCAPS



#subclass for tacs steady interface
class wedgeTACS(TacsSteadyInterface):
    def __init__(self, comm, tacs_comm, model, n_tacs_procs, datFile):
        super(wedgeTACS,self).__init__(comm, tacs_comm, model)

        assembler = None
        self.tacs_proc = False
        if comm.Get_rank() < n_tacs_procs:
            self.tacs_proc = True

            # Instantiate FEASolver
            structOptions = {
                'printtiming':True,
            }

            FEASolver = pyTACS(datFile, options=structOptions, comm=tacs_comm)

            # Material properties
            rho = 2780.0        # density kg/m^3
            E = 73.1e9          # Young's modulus (Pa)
            nu = 0.33           # Poisson's ratio
            kcorr = 5.0/6.0     # shear correction factor
            ys = 324.0e6        # yield stress
            specific_heat = 920.096
            cte = 24.0e-6
            kappa = 230.0

            tInput = 0.001*np.ones(3)

            def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
                elemIndex = kwargs['propID'] - 1
                t = tInput[elemIndex]
                prop = constitutive.MaterialProperties(rho=rho, specific_heat=specific_heat,
                                                       E=E, nu=nu, ys=ys, cte=cte, kappa=kappa)
                con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dvNum)

                elemList = []
                transform = None
                for elemDescript in elemDescripts:
                    if elemDescript in ['CQUAD4', 'CQUADR']:
                        elem = elements.Quad4ThermalShell(transform, con)
                    else:
                        print("Uh oh, '%s' not recognized" % (elemDescript))
                    elemList.append(elem)

                # Add scale for thickness dv
                scale = [1.0]
                return elemList, scale

            # Set up elements and TACS assembler
            FEASolver.initialize(elemCallBack)
            assembler = FEASolver.assembler

        self._initialize_variables(assembler, thermal_index=6)
        self.initialize(model.scenarios[0],model.bodies)

    def post_export_f5(self):
        flag = (TACS.OUTPUT_CONNECTIVITY |
                TACS.OUTPUT_NODES |
                TACS.OUTPUT_DISPLACEMENTS |
                TACS.OUTPUT_STRAINS |
                TACS.OUTPUT_STRESSES |
                TACS.OUTPUT_EXTRAS)
        f5 = TACS.ToFH5(self.assembler, TACS.BEAM_OR_SHELL_ELEMENT, flag)
        f5.writeToFile('TACSoutput.f5')

#class for NACA OML optimization, full aerothermoelastic, with ESP/CAPS parametric geometries
class NacaOMLOptimization():
    def __init__(self, comm, structCSM, fluidCSM, DVdict, analysisType = "aerothermoelastic"):

        self.comm = comm

        #status file
        if (self.comm.Get_rank() == 0):
            prevStatus = os.path.join(os.getcwd(), "prev_status.txt")
            statusFile = os.path.join(os.getcwd(), "status.txt")
            if (os.path.exists(statusFile)): os.remove(statusFile)
            self.status =  open(statusFile, "w")

        self.cwrite("Aerothermoelastic Optimization with FuntoFem and ESP/CAPS\n")
        self.cwrite("\tDesign Problem: NACA Symmetric Wing\n")
        self.cwrite("Authors: Sean Engelstad, Sejal Sahu, Graeme Kennedy\n")
        self.cwrite("\tGeorgia Tech SMDO Lab April 2022\n")
        self.cwrite("----------------------------")
        self.cwrite("----------------------------\n")


        #iteration counter for optimizer
        self.iteration = 1

        self.curDir = os.getcwd()

        #design variables dictionary
        #"name" : dvname
        #"type" : "thick", "shape", "aero"
        #"value" : initially zero
        self.DVdict = DVdict

        #type of analysis
        self.analysis_type = analysisType

        #load pointwise
        #module load pointwise/18.5R1
        #need to run module load pointwise/18.5R1 in terminal for it to work

        #initialize AIMS
        if (self.comm.Get_rank() == 0): self.initializeAIMs(structCSM, fluidCSM)

    def cwrite(self, text):
        if (self.comm.Get_rank() == 0):
            #write to the status file
            self.status.write(text)
        
            #immediately update it to be visibile in the file
            self.status.flush()

    def initializeAIMs(self, structCSM, fluidCSM):
        if (self.comm.Get_rank() == 0):
            #initialize all 6 ESP/CAPS AIMs used for the fluid and structural analysis

            #initialize pyCAPS structural problem
            self.capsStruct = pyCAPS.Problem(problemName = "struct",
                        capsFile = structCSM,
                        outLevel = 1)
            self.cwrite("Initialized caps Struct AIM\n")

            #initialize pyCAPS fluid problem
            self.capsFluid = pyCAPS.Problem(problemName = "fluid",
                        capsFile = fluidCSM,
                        outLevel = 1)
            self.cwrite("Initialized caps fluid AIM\n")

            #initialize egads Aim
            self.egadsAim = self.capsStruct.analysis.create(aim="egadsTessAIM")
            self.cwrite("Initialized egads AIM\n")

            #initialize tacs Aim
            self.tacsAim = self.capsStruct.analysis.create(aim = "tacsAIM", name = "tacs")
            self.cwrite("Initialized tacs AIM\n")

            #initialize pointwise AIM
            self.pointwiseAim = self.capsFluid.analysis.create(aim = "pointwiseAIM",
                                        name = "pointwise")
            self.cwrite("Initialized pointwise AIM\n")

            #initialize FUN3D AIM from Pointwise mesh
            self.fun3dAim = self.capsFluid.analysis.create(aim = "fun3dAIM",
                                    name = "fun3d")
            self.cwrite("Initialized fun3d AIM\n")

            #structural mesh settings
            self.structureMeshSettings()
            self.cwrite("Set Structure mesh settings\n")

            #fluid mesh settings
            self.fluidMeshSettings()
            self.cwrite("Set Fluid mesh settings\n")

            #set fun3d settings
            self.fun3dSettings()
            self.cwrite("Set fun3d settings\n")

    def structureMeshSettings(self):
        #Egads Aim section, for mesh
        self.egadsAim.input.Edge_Point_Min = 5
        self.egadsAim.input.Edge_Point_Max = 10

        self.egadsAim.input.Mesh_Elements = "Quad"

        self.egadsAim.input.Tess_Params = [.25,.01,15]

        #increase the precision in the BDF file
        self.tacsAim.input.File_Format = "Large"
        self.tacsAim.input.Mesh_File_Format = "Large"

        # Link the mesh
        self.tacsAim.input["Mesh"].link(self.egadsAim.output["Surface_Mesh"])

        # Set analysis type
        self.tacsAim.input.Analysis_Type = "Static"

        #materials section    
        madeupium    = {"materialType" : "isotropic",
                        "youngModulus" : 72.0E9 ,
                        "poissonRatio": 0.33,
                        "density" : 2.8E3,
                        "tensionAllow" :  20.0e7}

        self.tacsAim.input.Material = {"madeupium": madeupium}

        # Material properties section
        OMLshell = {"propertyType" : "Shell",
                    "membraneThickness" : 0.01,
                    "material"        : "madeupium",
                    "bendingInertiaRatio" : 1.0, # Default
                    "shearMembraneRatio"  : 5.0/6.0} # Default
        ribshell  = {"propertyType" : "Shell",
                    "membraneThickness" : 0.02,
                    "material"        : "madeupium",
                    "bendingInertiaRatio" : 1.0, # Default
                    "shearMembraneRatio"  : 5.0/6.0} # Default

        sparshell = {"propertyType" : "Shell",
                    "membraneThickness" : 0.05,
                    "material"        : "madeupium",
                    "bendingInertiaRatio" : 1.0, # Default
                    "shearMembraneRatio"  : 5.0/6.0} # Default

        self.tacsAim.input.Property = {"rib": ribshell,
        "spar" : sparshell,
        "OML" : OMLshell}

        # constraint section
        constraint1 = {"groupName" : "wingRoot",
                        "dofConstraint" : 123456}

        self.tacsAim.input.Constraint = {"fixRoot": constraint1}

        # don't set any loads

        ##-----------setup DVs and DVRs-----------##

        #list of design variables, with thick1-thickN for thickness DVs
        desvars = ["area","aspect","taper","ctwist","lesweep","dihedral","thick1", "thick2", "thick3"]
        nvar = len(desvars)

        #where the thickDVs have thick1 is for "rib", thick2 for "spar" etc.
        thickness = 0.01

        #make initial DV and DVR dict
        DVdict = {}
        DVRdict = {}

        #add thickDVs and geomDVs to caps
        thickCt = 0
        for DV in self.DVdict:
            dvname = DV["name"]
            if (DV["type"] == "thick"):
                #add thickDV entry into DV_Relations and DV Dicts
                DVRdict[dvname] = self.makeThicknessDVR(dvname)
                DVdict[dvname] = self.makeThicknessDV(DV["capsGroup"],thickness)
                thickCt += 1
            elif (DV["type"] == "shape"): #geomDV, add empty entry into DV dicts
                DVdict[dvname] = {}

        #input DVdict and DVRdict into tacsAim
        self.tacsAim.input.Design_Variable = DVdict
        self.tacsAim.input.Design_Variable_Relation = DVRdict

    def fluidMeshSettings(self):
        # Dump VTK files for visualization
        self.pointwiseAim.input.Proj_Name   = "TransportWing"
        self.pointwiseAim.input.Mesh_Format = "VTK"

        # Connector level
        self.pointwiseAim.input.Connector_Turn_Angle       = 10
        self.pointwiseAim.input.Connector_Prox_Growth_Rate = 1.2
        self.pointwiseAim.input.Connector_Source_Spacing   = True

        # Domain level
        self.pointwiseAim.input.Domain_Algorithm    = "AdvancingFront"
        self.pointwiseAim.input.Domain_Max_Layers   = 15
        self.pointwiseAim.input.Domain_Growth_Rate  = 1.25
        self.pointwiseAim.input.Domain_TRex_ARLimit = 40.0
        self.pointwiseAim.input.Domain_Decay        = 0.8
        self.pointwiseAim.input.Domain_Iso_Type = "Triangle"

        # Block level
        self.pointwiseAim.input.Block_Boundary_Decay       = 0.8
        self.pointwiseAim.input.Block_Collision_Buffer     = 1.0
        self.pointwiseAim.input.Block_Max_Skew_Angle       = 160.0
        self.pointwiseAim.input.Block_Edge_Max_Growth_Rate = 1.5
        self.pointwiseAim.input.Block_Full_Layers          = 1
        self.pointwiseAim.input.Block_Max_Layers           = 100
        self.pointwiseAim.input.Block_TRexType = "TetPyramid"
        #T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms).

        # Set wall spacing for capsMesh == leftWing and capsMesh == riteWing
        viscousWall  = {"boundaryLayerSpacing" : 0.001}
        self.pointwiseAim.input.Mesh_Sizing = {"OML": viscousWall,
                "Farfield": {"bcType":"Farfield"},
                "Symmetry" : {"bcType" : "SymmetryY"}}


    def fun3dSettings(self):
        self.fun3dAim.input.Boundary_Condition = {"OML": {"bcType" : "Inviscid"},
                "Farfield": {"bcType":"Farfield"},
                "Symmetry" : "SymmetryY"}

        # # Set project name
        # fun3dAIM.input.Proj_Name = "fun3dTetgenTest"
        
        # # Link the mesh
        # fun3dAIM.input["Mesh"].link(meshAIM.output["Volume_Mesh"])
        
        # fun3dAIM.input.Mesh_ASCII_Flag = False
        
        # # Set AoA number
        # myProblem.analysis["fun3d"].input.Alpha = 1.0
        
        # # Set Mach number
        # myProblem.analysis["fun3d"].input.Mach = 0.5901
        
        # # Set equation type
        # fun3dAIM.input.Equation_Type = "compressible"
        
        # # Set Viscous term
        # myProblem.analysis["fun3d"].input.Viscous = "inviscid"
        
        # # Set number of iterations
        # myProblem.analysis["fun3d"].input.Num_Iter = 10
        
        # # Set CFL number schedule
        # myProblem.analysis["fun3d"].input.CFL_Schedule = [0.5, 3.0]
        
        # # Set read restart option
        # fun3dAIM.input.Restart_Read = "off"
        
        # # Set CFL number iteration schedule
        # myProblem.analysis["fun3d"].input.CFL_Schedule_Iter = [1, 40]
        
        # # Set overwrite fun3d.nml if not linking to Python library
        # myProblem.analysis["fun3d"].input.Overwrite_NML = True

    def initF2F(self, x):
        structDVs = x["struct"]
        
        maximum_mass = 40.0 
        num_tacs_dvs = len(structDVs)

        # Set up the communicators
        n_tacs_procs = 1

        world_rank = self.comm.Get_rank()
        if world_rank < n_tacs_procs:
            color = 55
            key = world_rank
        else:
            color = MPI.UNDEFINED
            key = world_rank
        tacs_comm = self.comm.Split(color,key)

        #==================================================================================================#
        # Originally _build_model()

        # Build the model
        self.model = FUNtoFEMmodel('NACA Wing Simulation')
        self.wing = Body('wing', analysis_type=self.analysis_type, group=0,boundary=1)

        for i in range(num_tacs_dvs):
            self.wing.add_variable('structural',Variable('thickness '+ str(i),value=structDVs[i],lower = 0.0001, upper = 0.01))

        self.model.add_body(self.wing)

        steady = Scenario('steady', group=0, steps=5)
        function1 = Function('ksfailure',analysis_type='structural')
        steady.add_function(function1)

        function2 = Function('mass',analysis_type='structural', adjoint=False) #,adjoint=False
        steady.add_function(function2)

        # function3 = Function('lift', analysis_type='aerodynamic')
        # steady.add_function(function3)

        # function4 = Function('drag', analysis_type='aerodynamic')
        # steady.add_function(function4)

        self.model.add_scenario(steady)

        self.cwrite("setup model, ")

        #==================================================================================================#

        # instantiate TACS on the master
        solvers = {}
        solvers['flow'] = Fun3dInterface(self.comm,self.model,flow_dt=1.0)
        self.cwrite("setup fun3d interface, ")
        datFile = os.path.join(self.curDir,"steady","Flow","nastran_CAPS.dat")
        solvers['structural'] = wedgeTACS(self.comm,tacs_comm,self.model,n_tacs_procs, datFile)
        self.cwrite("setup tacs interface\n")

        # L&D transfer options
        transfer_options = {'analysis_type': self.analysis_type,
                            'scheme': 'meld', 'thermal_scheme': 'meld'}

        # instantiate the driver
        self.driver = FUNtoFEMnlbgs(solvers,self.comm,tacs_comm,0,self.comm,0,transfer_options,model=self.model)
        struct_tacs = solvers['structural'].assembler
        self.cwrite("\t setup adjoint driver, ")

        #section to update funtofem DV values for next fun3d run
        self.model.set_variables(structDVs)

    def forwardAnalysis(self, x):
        #update status
        self.cwrite("----------------------------")
        self.cwrite("----------------------------\n")
        self.cwrite("Iteration #{}\n".format(self.iteration))
        self.iteration += 1

        #set design variables from desvarDict
        self.updateDesign(x)
        thickDVstr = "[rib,spar,OML]"
        shapeDVstr = "[area,aspect,camb0,cambf,ctwist,dihedral,lesweep,taper,tc0,tcf]"
        self.cwrite("thickDVs {}\n\t{}\nshapeDVs {}\n\t{}\n".format(thickDVstr, x["struct"],shapeDVstr, x["shape"]))

        #generate structure mesh with egads and tacs AIMs
        self.buildStructureMesh()
        self.cwrite("built structure mesh\n")

        #generate fluid mesh with Pointwise
        self.buildFluidMesh()
        self.cwrite("built fluid mesh\n")

        #initialize fun2fem with new meshes
        self.cwrite("Initializing funtofem... ")
        self.initF2F(x)
        self.cwrite("initialized funtofem\n")

        #run FUNtoFEM forward analysis
        #run funtofem forward analysis, including fun3d
        self.cwrite("Running F2F forward analysis... ")
        self.driver.solve_forward()

        #get function values from forward analysis
        functions = self.model.get_functions()

        self.cwrite("completed F2F forward analysis\n")

        #send or store function values
        return functions

    def adjointAnalysis(self, x):
        #update status
        self.cwrite("Running F2F adjoint analysis... ")

        #run funtofem adjoint analysis, with fun3d
        self.driver.solve_adjoint()

        #update status
        self.cwrite("finished adjoint analysis\n")

        #get the funtofem gradients
        f2fgrads = self.model.get_function_gradients()
        
        #update gradient for non shape DVs
        structGrad = np.zeros(3)
        ct = 0
        for funcKey in self.funcKeys:
            for DV in self.DVdict:
                dvname = DV["name"]
                if (not(DV["type"] == "shape")): structGrad[ct] = f2fgrads[funcKey][ct]
                ct += 1
        
        #get aero and struct mesh sensitivities for shape DVs
        self.aeroIds, self.aero_mesh_sens = self.wing.collect_coordinate_derivatives(self.comm, "aero")
        self.structIds, self.struct_mesh_sens = self.wing.collect_coordinate_derivatives(self.comm, "struct")

        #compute shape derivatives from aero and struct mesh sensitivities
        self.computeShapeDerivatives()
        self.cwrite("computed shape derivative chain rule products\n")


        #send or store gradient
        return structGrad, shapeGrad

    def updateDesign(self, x):

        thickDVs = x["struct"]
        shapeDVs = x["shape"]

        #update the design variable dictionaries with new DV values
        thickCt = 0
        shapeCt = 0
        ct = 0
        for DVdict in self.DVdict:
            if (DVdict["type"] == "thick"):
                DVdict["value"] = thickDVs[thickCt]
                self.DVdict[ct] = DVdict
                thickCt += 1
                ct += 1
            elif (DVdict["type"] == "shape"):
                DVdict["value"] = shapeDVs[shapeCt]
                self.DVdict[ct] = DVdict
                shapeCt += 1
                ct += 1

        
        #update shapeDVs in each caps problem
        #update thickness design variables in tacsAim
        if (self.comm.Get_rank() == 0):
            #grab the dictionaries in tacsAim
            propDict = self.tacsAim.input.Property
            DVRdict = self.tacsAim.input.Design_Variable_Relation
            DVdict = self.tacsAim.input.Design_Variable

            #update thickness design variables in caps aims
            thickCt = 0
            thickVec = []
            for DV in self.DVdict:
                dvname = DV["name"]
                value = DV["value"]
                if (DV["type"] == "shape"):
                    #update both of the capsStruct and capsFluid geometries
                    self.capsStruct.geometry.despmtr[dvname].value = value
                    self.capsFluid.geometry.despmtr[dvname].value = value
                elif (DV["type"] == "thick"):
                    #update the thickness DV in its tacsAim dictionaries
                    DVRdict[dvname] = self.makeThicknessDVR(dvname)
                    capsGroup = DVdict[dvname]["groupName"]
                    DVdict[dvname] = self.makeThicknessDV(capsGroup,value)
                    propDict[capsGroup]["membraneThickness"] = value

                    #update funtofem thickness vec
                    thickVec.append(value)

            #update tacsAim dictionaries
            self.tacsAim.input.Property = propDict
            self.tacsAim.input.Design_Variable_Relation = DVRdict
            self.tacsAim.input.Design_Variable = DVdict
        
    def makeThicknessDV(self, capsGroup, thickness):
        #thick DV dictionary for Design_Variable Dict
        desvar    = {"groupName" : capsGroup,
                "initialValue" : thickness,
                "lowerBound" : thickness*0.5,
                "upperBound" : thickness*1.5,
                "maxDelta"   : thickness*0.1}
        return desvar
        
    def makeThicknessDVR(self, DVname):
        #thick DV dictionary for Design_Variable_Relation Dict
        DVR = {"variableType": "Property",
        "fieldName" : "T",
        "constantCoeff" : 0.0,
        "groupName" : DVname,
        "linearCoeff" : 1.0}
        return DVR           

    def buildStructureMesh(self):
        self.cwrite("Building structure mesh... ")
        #build structure mesh by running tacsAim preanalysis
        if (self.comm.Get_rank() == 0):
            self.tacsAim.preAnalysis()

            #move struct mesh files to flow directory 
            for extension in [".bdf",".dat"]:
                src = os.path.join(self.tacsAim.analysisDir, self.tacsAim.input.Proj_Name+extension)
                dest = os.path.join(self.curDir,"steady","Flow",self.tacsAim.input.Proj_Name+extension)
                shutil.copy(src, dest)

    def buildFluidMesh(self):
        self.cwrite("Building fluid mesh... ")

        if (self.comm.Get_rank() == 0):
            #build fluid mesh by running pointwise and then linking with fun3d
            self.runPointwise()
            self.cwrite("ran pointwise, ")


            #update fun3d with current mesh so it knows mesh sensitivity
            self.fun3dAim.input["Mesh"].link(self.pointwiseAim.output["Volume_Mesh"])
            self.fun3dAim.preAnalysis()
            self.cwrite("linked to fun3d, ")

            #copy the mesh file to funtofem directory steady/flow/
            filename = "caps.GeomToMesh.ugrid"
            src = os.path.join(self.pointwiseAim.analysisDir, "caps.GeomToMesh.ugrid")
            dest = os.path.join(self.curDir, "steady", "Flow", "caps.GeomToMesh.ugrid")
            shutil.copy(src, dest)

    def runPointwise(self):
        #run AIM pre-analysis
        self.pointwiseAim.preAnalysis()

        #move to test directory
        self.curDir = os.getcwd()
        os.chdir(self.pointwiseAim.analysisDir)

        CAPS_GLYPH = os.environ["CAPS_GLYPH"]
        for i in range(1): #can run extra times if having license issues
                os.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")
            #if os.path.isfile('caps.GeomToMesh.gma') and os.path.isfile('caps.GeomToMesh.ugrid'): break

        os.chdir(self.curDir)

        #run AIM postanalysis, files in self.pointwiseAim.analysisDir
        self.pointwiseAim.postAnalysis()      


    def objCon(self, x):
        functions = self.forwardAnalysis(x)
        #0 - stress
        #1 - mass
        #2 - lift
        #3 - drag

        funcs = {}
        funcs["obj"] = functions[1] #lift
        #funcs["con"] = func2
        fail = False

        #update status
        self.cwrite("Objective function is {}\n".format(funcs["obj"]))

        return funcs, fail

    def objGrad(self, x):
        gradients = self.adjointAnalysis(x)
        #0 - stress
        #1 - mass
        #2 - lift
        #3 - drag
        fail = False
        grad1 = gradients[1][:]

        sens = {}
        sens["obj"] = { "struct": structGrad,
                        "shape" : shapeGrad}
        #sens["con"] = {"x": [grad2]}

        #update status
        self.cwrite("gradient is {}\n".format(sens["obj"]))
        return objGrad, fail

    def computeShapeDerivatives(self):
        if (self.comm.Get_rank() == 0):
            #add struct_mesh_sens part to shape DV derivatives#struct shape derivatives
            self.applyStructMeshSens()

            #add aero_mesh_sens part to shape DV derivatives
            self.applyAeroMeshSens()

    def applyStructMeshSens(self):
        #print struct mesh sens to struct mesh sens file
        #where to print .sens file
        structSensFile = os.path.join(self.tacsAim.analysisDir, self.tacsAim.input.Proj_Name+".sens")
        
        #open the file
        with open(structSensFile, "w") as f:
            
            #write (nfunctions) in first line
            f.write("{}\n".format(self.nfunc))
            
            funcInd = 0
            #for each function mass, stress, etc.
            for key in self.funcKeys:
                
                #get the pytacs/tacs sensitivity w.r.t. mesh for that function
                sens = self.struct_mesh_sens[:,funcInd]

                #write the key,value,nnodes of the function
                f.write(key + "\n")
                f.write("{}\n".format(self.function[key]))
                f.write("{}\n".format(self.structNodes))

                #for each node, print nodeind, dfdx, dfdy, dfdz for that mesh element
                for structId in self.structIds: # d(Func1)/d(xyz)
                    bdfind = structId
                    #bdfind = nodeind + 1
                    f.write("{} {} {} {}\n".format(bdfind, sens[nodeind,0], sens[nodeind,1], sens[nodeind,2]))

                funcInd += 1
        
        #update status
        self.cwrite("printed struct.sens file\n")

        #run aim postanalysis
        self.tacsAim.postAnalysis()
        self.cwrite("completed tacsAim postAnalysis()\n")

        #update shape DV derivatives from struct mesh part
        for funcKey in self.funcKeys:
            for DV in self.DVdict:
                dvname = DV["name"]
                if (DV["type"] == "shape"): self.gradient[funcKey][dvname] = self.tacsAim.dynout[funcKey].deriv(dvname)

        #update status
        self.cwrite("finished shape DV contribution from struct mesh sens\n")


    def applyAeroMeshSens(self):
        #print aero mesh sens to aero mesh sens file
        #where to print .sens file
        structSensFile = os.path.join(self.fun3dAim.analysisDir, self.fun3dAim.input.Proj_Name+".sens")
        
        #open the file
        with open(structSensFile, "w") as f:
            
            #write (nfunctions) in first line
            f.write("{}\n".format(self.nfunc))
            
            funcInd = 0
            #for each function mass, stress, etc.
            for key in self.funcKeys:
                
                #get the pytacs/tacs sensitivity w.r.t. mesh for that function
                sens = self.aero_mesh_sens[:,funcInd]

                #write the key,value,nnodes of the function
                f.write(key + "\n")
                f.write("{}\n".format(self.function[key]))
                f.write("{}\n".format(self.aeroNodes))

                #for each node, print nodeind, dfdx, dfdy, dfdz for that mesh element
                for aeroId in self.aeroIds: # d(Func1)/d(xyz)
                    bdfind = aeroId
                    #bdfind = nodeind + 1
                    f.write("{} {} {} {}\n".format(bdfind, sens[nodeind,0], sens[nodeind,1], sens[nodeind,2]))

                funcInd += 1
        
        #update status
        self.cwrite("printed aero.sens file\n")

        #run aim postanalysis
        self.fun3dAim.postAnalysis()
        self.cwrite("completed pointwiseAim postAnalysis()\n")

        #update shape DV derivatives from aero mesh part
        for funcKey in self.funcKeys:
            for DV in self.DVdict:
                dvname = DV["name"]
                if (DV["type"] == "shape"): self.gradient[funcKey][dvname] += self.fun3dAim.dynout[funcKey].deriv(dvname)
            
        #update status
        self.cwrite("finished shape DV contribution from aero mesh sens\n")


##----------Outside of class, run cases------------------##

#setup the design variables
#make shape DVs
DVdict = []
for dvname in ["area","aspect","camb0","cambf","ctwist", "dihedral","lesweep", "taper","tc0","tcf"]:
    tempDict = {"name" : dvname,
                "type" : "shape",
                "value" : 0.0,
                "capsGroup" : ""}
    DVdict.append(tempDict)

#setup thick DVs
capsGroup = ["rib","spar","OML"]
thickCt = 0
for dvname in ["thick1", "thick2", "thick3"]:
    tempDict = {"name" : dvname,
                "type" : "thick",
                "value" : 0.1,
                "capsGroup" : capsGroup[thickCt]}
    thickCt += 1
    DVdict.append(tempDict)

#call the class and initialize it
comm = MPI.COMM_WORLD

nacaOpt = NacaOMLOptimization(comm, "naca_OML_struct.csm", "naca_OML_fluid.csm", DVdict, "aerothermoelastic")

#setup pyOptSparse
sparseProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", nacaOpt.objCon)

#  ["area","aspect","camb0","cambf","ctwist", "dihedral","lesweep", "taper","tc0","tcf"]
lbnds =   [20.0, 3.0, 0.0 , 0.0, 1.0,  1.0, 3.0,   0.3, 0.0, 0.0]
init = [40.0, 6.0,  0.0, 0.0, 5.0,  5.0, 30.0,  0.5, 0.1, 0.1]
ubnds =   [100.0,10.0, 0.3, 0.3, 10.0, 20.0, 50.0, 1.0, 0.3, 0.3]

#["rib","spar","OML"]
lBnds2 = 0.0001 * np.ones(3)
uBnds2 = 0.01*np.ones(3)
init2 = 0.001*np.ones(3)

sparseProb.addVarGroup("shape", 10, "c", lower=lbnds, upper=ubnds, value=init)
sparseProb.addVarGroup("struct", 3, "c", lower=lBnds2, upper=uBnds2, value=init2)

#optProb.addConGroup("con", 1, lower=1, upper=1)
sparseProb.addObj("obj")

# if comm.rank == 0:
#     print(sparseProb)

optOptions = {"IPRINT": -1}
opt = SLSQP(options=optOptions)
sol = opt(sparseProb, sens=nacaOpt.objGrad)

if comm.rank == 0:
    print(sol)
    print('\nsol.xStar:  ', sol.xStar)

#close the status file
nacaOpt.status.close()