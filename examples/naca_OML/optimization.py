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
import os
import numpy as np

#import other modules
from pyoptsparse import SLSQP, Optimization
import pyCAPS



#subclass for tacs steady interface
class wedgeTACS(TacsSteadyInterface):
    def __init__(self, comm, tacs_comm, model, n_tacs_procs):
        super(wedgeTACS,self).__init__(comm, tacs_comm, model)

        assembler = None
        self.tacs_proc = False
        if comm.Get_rank() < n_tacs_procs:
            self.tacs_proc = True

            # Instantiate FEASolver
            structOptions = {
                'printtiming':True,
            }

            bdfFile = os.path.join(os.path.dirname(__file__), 'nastran_CAPS.dat')
            FEASolver = pyTACS(bdfFile, options=structOptions, comm=tacs_comm)

            # Material properties
            rho = 2780.0        # density kg/m^3
            E = 73.1e9          # Young's modulus (Pa)
            nu = 0.33           # Poisson's ratio
            kcorr = 5.0/6.0     # shear correction factor
            ys = 324.0e6        # yield stress
            specific_heat = 920.096
            cte = 24.0e-6
            kappa = 230.0

            t = 0.001

            def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
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
    def __init__(self, structCSM, fluidCSM, DVdict, analysisType = "aerothermoelastic"):

        #design variables dictionary
        #"name" : dvname
        #"type" : "thick", "shape", "aero"
        #"capsGroup" : corr to thick, otherwise empty
        #"value" : initially zero
        self.DVdict = DVdict
        self.analysis_type = analysisType

        #load pointwise
        #need to run module load pointwise/18.5R1 in terminal for it to work

        #initialize AIMS
        self.initializeAIMs(structCSM, fluidCSM)

        #initialize funtofem
        self.initFun2Fem()

    def initializeAIMs(self, structCSM, fluidCSM):
        #initialize all 6 ESP/CAPS AIMs used for the fluid and structural analysis

        #initialize pyCAPS structural problem
        self.capsStruct = pyCAPS.Problem(problemName = "struct",
                    capsFile = structCSM,
                    outLevel = 1)

        #initialize pyCAPS fluid problem
        self.capsFluid = pyCAPS.Problem(problemName = "fluid",
                    capsFile = fluidCSM,
                    outLevel = 1)

        #initialize egads Aim
        self.egadsAim = self.capsStruct.analysis.create(aim="egadsTessAIM")

        #initialize tacs Aim
        self.tacsAim = self.capsStruct.analysis.create(aim = "tacsAIM", name = "tacs")  

        #initialize pointwise AIM
        self.pointwiseAIM = self.capsFluid.analysis.create(aim = "pointwiseAIM",
                                    name = "pointwise")

        #initialize FUN3D AIM from Pointwise mesh
        self.fun3dAIM = self.capsFluid.analysis.create(aim = "fun3dAIM",
                                name = "fun3d")

        #structural mesh settings
        self.structureMeshSettings()

        #fluid mesh settings
        self.fluidMeshSettings()

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
        self.pointwiseAIM.input.Proj_Name   = "TransportWing"
        self.pointwiseAIM.input.Mesh_Format = "VTK"

        # Connector level
        self.pointwiseAIM.input.Connector_Turn_Angle       = 10
        self.pointwiseAIM.input.Connector_Prox_Growth_Rate = 1.2
        self.pointwiseAIM.input.Connector_Source_Spacing   = True

        # Domain level
        self.pointwiseAIM.input.Domain_Algorithm    = "AdvancingFront"
        self.pointwiseAIM.input.Domain_Max_Layers   = 15
        self.pointwiseAIM.input.Domain_Growth_Rate  = 1.25
        self.pointwiseAIM.input.Domain_TRex_ARLimit = 40.0
        self.pointwiseAIM.input.Domain_Decay        = 0.8
        self.pointwiseAIM.input.Domain_Iso_Type = "Triangle"

        # Block level
        self.pointwiseAIM.input.Block_Boundary_Decay       = 0.8
        self.pointwiseAIM.input.Block_Collision_Buffer     = 1.0
        self.pointwiseAIM.input.Block_Max_Skew_Angle       = 160.0
        self.pointwiseAIM.input.Block_Edge_Max_Growth_Rate = 1.5
        self.pointwiseAIM.input.Block_Full_Layers          = 1
        self.pointwiseAIM.input.Block_Max_Layers           = 100
        self.pointwiseAIM.input.Block_TRexType = "TetPyramid"
        #T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms).

        # Set wall spacing for capsMesh == leftWing and capsMesh == riteWing
        viscousWall  = {"boundaryLayerSpacing" : 0.001}
        self.pointwiseAIM.input.Mesh_Sizing = {"OML": viscousWall,
                "Farfield": {"bcType":"Farfield"},
                "Symmetry" : {"bcType" : "SymmetryY"}}

    def initFun2Fem(self):
        analysis_type='aerothermoelastic'

        maximum_mass = 40.0 
        num_tacs_dvs = 3

        # Set up the communicators
        n_tacs_procs = 1

        comm = MPI.COMM_WORLD
        self.comm = comm

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

        thickness = 0.001

        # Build the model
        self.model = FUNtoFEMmodel('NACA Wing Simulation')
        wing = Body('wing', analysis_type=analysis_type, group=0,boundary=1)

        for i in range(num_tacs_dvs):
            wing.add_variable('structural',Variable('thickness '+ str(i),value=thickness,lower = 0.0001, upper = 0.01))

        self.model.add_body(wing)

        steady = Scenario('steady', group=0, steps=5)
        function1 = Function('ksfailure',analysis_type='structural')
        steady.add_function(function1)

        function2 = Function('mass',analysis_type='structural',adjoint=False)
        steady.add_function(function2)

        function3 = Function('lift', analysis_type='aerodynamic')
        steady.add_function(function3)

        function4 = Function('drag', analysis_type='aerodynamic')
        steady.add_function(function4)
        
        self.model.add_scenario(steady)

        #==================================================================================================#

        # instantiate TACS on the master
        solvers = {}
        solvers['flow'] = Fun3dInterface(comm,model,flow_dt=1.0)
        solvers['structural'] = wedgeTACS(comm,tacs_comm,self.model,n_tacs_procs)

        # L&D transfer options
        transfer_options = {'analysis_type': analysis_type,
                            'scheme': 'meld', 'thermal_scheme': 'meld'}

        # instantiate the driver
        self.driver = FUNtoFEMnlbgs(solvers,self.comm,tacs_comm,0,comm,0,transfer_options,model=self.model)
        struct_tacs = solvers['structural'].assembler

        obj_scale = 0.0106
        con_scale = 3.35

    def forwardAnalysis(self, x):
        #set design variables from desvarDict
        self.updateDesign(x)

        #generate structure mesh with egads and tacs AIMs
        self.buildStructureMesh()

        #generate fluid mesh with Pointwise
        self.buildFluidMesh()

        #run FUNtoFEM forward analysis
        #remember forward solution
        self.fun2femForward()

        #send or store function values
        return self.function

    def adjointAnalysis(self, x):
        #run FUNtoFEM adjoint analysis
        self.fun2femAdjoint()

        #compute shape derivatives from aero and struct mesh sensitivities
        self.computeShapeDerivatives()

        #send or store gradient
        return self.gradient

    def updateDesign(self, desvars):
        #update the design variable dictionaries with new DV values
        ct = 0
        for DVdict in self.DVdict:
            DVdict["value"] = desvars[ct]
            self.DVdict[ct] = DVdict
            ct += 1
        
        #update shapeDVs in each caps problem
        #update thickness design variables in tacsAim

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
                DVdict[deskey] = self.makeThicknessDV(capsGroup,value)
                propDict[capsGroup]["membraneThickness"] = v

                #update funtofem thickness vec
                thickVec.append(value)

        #update tacsAim dictionaries
        self.tacsAim.input.Property = propDict
        self.tacsAim.input.Design_Variable_Relation = DVRdict
        self.tacsAim.input.Design_Variable = DVdict

        #section to update funtofem aero DV values for next fun3d run
        self.model.set_variables(np.array(thickVec))
        
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
        #build structure mesh by running tacsAim preanalysis
        self.tacsAim.preAnalysis()

    def buildFluidMesh(self):
        #build fluid mesh by running pointwise and then linking with fun3d
        self.runPointwise()

        #update fun3d with current mesh so it knows mesh sensitivity
        self.fun3dAIM.input["Mesh"].link(self.pointwiseAim.output["Volume_Mesh"])
        self.fun3dAim.preAnalysis()

    def runPointwise(self):
        #run AIM pre-analysis
        self.pointwiseAim.preAnalysis()

        #move to test directory
        currentDir = os.getcwd()
        os.chdir(self.pointwiseAim.analysisDir)

        CAPS_GLYPH = os.environ["CAPS_GLYPH"]
        for i in range(1): #can run extra times if having license issues
                os.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")
            #if os.path.isfile('caps.GeomToMesh.gma') and os.path.isfile('caps.GeomToMesh.ugrid'): break

        os.chdir(currentDir)

        #run AIM postanalysis, files in self.pointwiseAim.analysisDir
        self.pointwiseAim.postAnalysis()

    def fun2femForward(self):        
        #run funtofem forward analysis, including fun3d
        self.driver.solve_forward()

        #get function values from forward analysis
        self.functions = self.model.get_functions()

    def fun2femAdjoint(self):
        #run funtofem adjoint analysis, with fun3d
        self.driver.solve_adjoint()

        #get the funtofem gradients
        f2fgrads = self.model.get_function_gradients()
        
        #update gradient for non shape DVs
        ct = 0
        for funcKey in self.funcKeys:
            for DV in self.DVdict:
                dvname = DV["name"]
                if (not(DV["type"] == "shape")): self.gradient[func][dvname] = f2fgrads[funcKey][ct]
                ct += 1
        

        #get aero and struct mesh sensitivities for shape DVs
        self.aero_mesh_sens = self.funtofem.body.aero_shape_term
        self.struct_mesh_sens = self.funtofem.body.struct_shape_term

    def objFunc(self, x):
        functions = self.forwardAnalysis(x)
        fail = 0

        return functions, fail

    def objGrad(self, x):
        gradients = self.adjointAnalysis(x)
        fail = 0
        return gradients, fail

    def computeShapeDerivatives(self):
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
                for nodeind in range(self.structNodes): # d(Func1)/d(xyz)
                    bdfind = nodeind + 1
                    f.write("{} {} {} {}\n".format(bdfind, sens[nodeind,0], sens[nodeind,1], sens[nodeind,2]))

                funcInd += 1
        
        #run aim postanalysis
        self.tacsAim.postAnalysis()
        print("ran postAnalysis()")

        #update shape DV derivatives from struct mesh part
        for funcKey in self.funcKeys:
            for DV in self.DVdict:
                dvname = DV["name"]
                if (DV["type"] == "shape"): self.gradient[dvname] = self.tacsAim.dynout[funcKey].deriv(dvname)


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
                for nodeind in range(self.aeroNodes): # d(Func1)/d(xyz)
                    bdfind = nodeind + 1
                    f.write("{} {} {} {}\n".format(bdfind, sens[nodeind,0], sens[nodeind,1], sens[nodeind,2]))

                funcInd += 1
        
        #run aim postanalysis
        self.fun3dAim.postAnalysis()
        print("ran postAnalysis()")

        #update shape DV derivatives from aero mesh part
        for funcKey in self.funcKeys:
            for DV in self.DVdict:
                dvname = DV["name"]
                if (DV["type"] == "shape"): self.gradient[dvname] += self.fun3dAim.dynout[funcKey].deriv(dvname)
            

##----------Outside of class, run cases------------------##

#setup the design variables
#make shape DVs
DVdict = []
for dvname in ["area","aspect","ctwist", "dihedral","lesweep", "taper"]:
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
nacaOpt = NacaOMLOptimization("naca_OML_struct.csm", "naca_OML_fluid.csm", DVdict, "aerothermoelastic")

#setup pyOptSparse
sparseProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", nacaOpt.objFunc)

sparseProb.addVarGroup("x", 12, "c", lower=0.0001*np.ones(12), upper=0.01*np.ones(12), value=0.001)

#optProb.addConGroup("con", 1, lower=1, upper=1)
optProb.addObj("obj")

comm = MPI.COMM_WORLD
# if comm.rank == 0:
print(sparseProb)

optOptions = {"IPRINT": -1}
opt = SLSQP(options=optOptions)
sol = opt(sparseProb, sens=nacaOpt.objGrad)

if comm.rank == 0:
    print(sol)
    print(sol, file=nacaOpt.optHist)
    print('\nsol.xStar:  ', sol.xStar)
    print('\nsol.xStar:  ', sol.xStar, file=sparseProb.optHist)
print(sol, file=sparseProb.optHistAll)
print('\nsol.xStar:  ', sol.xStar, file=sparseProb.optHistAll)

dp.optHist.close()
dp.optHistAll.close()