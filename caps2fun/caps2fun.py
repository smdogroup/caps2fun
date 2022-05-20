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


script to run funtofem from design variable inputs, funtofem.input file and .csm files
returns function, gradient information in funtofem.out file

Aerothermoelastic Optimization with FuntoFem and ESP/CAPS
	Design Problem: NACA Symmetric Wing
Authors: Sean Engelstad, Sejal Sahu, Graeme Kennedy
	Georgia Tech SMDO Lab April 2022
"""

#--------------------------- File Checklist --------------------------------------#

# Files to make sure you have before running F2F
# 1) fun3d.nml in steady/Flow folder - step# matches F2F, projname = meshStyle (pointwise/tetgen)
# 2) pointwise.mapbc if running pointwise mesher
# 3) a fluid CSM file for ESP/CAPS, fluid volume mesh
# 4) a struct CSM file for ESP/CAPS, the struct surf mesh

#check the DVdict was read in correctly
#print(DVdict)

#FUNtoFEM software is used here for full aerothermoelastic optimization of a NACA wing
#   with parametric geometries built by ESP/CAPS
#Authors: Sean Engelstad, Sejal Sahu, Graeme Kennedy
#Work for Aviation Paper 2022 on Full aerothermoelastic optimization with Parametric geometries of ESP/CAPS

from __future__ import print_function

#import normal classes
import os, shutil, sys
import time
import numpy as np
import pyCAPS
import f90nml

#read whether to use complex mode or not from input file
f2fin = os.path.join(os.getcwd(), "funtofem","funtofem.in")
isComplex = False
inputExists = os.path.exists(f2fin)
if (inputExists):
    hdl = open(f2fin,"r")
    lines = hdl.readlines()

    for line in lines:
        chunks = line.split(",")
        if ("mode" in line):
            isComplex = "complex_step" in line
    hdl.close()

#turn on complex mode if sent in through input file
#need to do this before importing pyfuntofem otherwise fun3d will import real flow solvers instead
if (isComplex): 
    os.environ['CMPLX_MODE'] = "1"
else: #otherwise set it to empty which means off
    os.environ["CMPLX_MODE"] = ""

#import funtofem classes and functions
from pyfuntofem.model  import *
from pyfuntofem.driver import *
from pyfuntofem.fun3d_interface import *
from pyfuntofem.tacs_interface import TacsSteadyInterface

#import from tacs
from tacs.pytacs import pyTACS
from tacs import functions
from tacs import TACS, functions, constitutive, elements, pyTACS, problems

#class to run Funtofem with full aerothermoelastic optimization, with parametric geometries from Engineering Sketch Pad
class Caps2Fun():
    def __init__(self):

        #initialize MPI
        comm = MPI.COMM_WORLD
        self.comm = comm

        #set root directory
        self.root_dir = os.getcwd()

        #read the config file
        self.readConfig()

        #read input from file
        self.readInput()

        #make the status file
        self.makeStatusFile()

        #load pointwise
        #module load pointwise/18.5R1
        #need to run module load pointwise/18.5R1 in terminal for it to work

        #initialize AIMS
        if (self.comm.Get_rank() == 0): self.initializeAIMs()

    def readConfig(self):
        #read the config files

        #initial values of each setting
        self.config = None

        #subfunction to parse the config file, used for each case
        def parseFile(file):
            handle = open(file, "r")
            lines = handle.readlines()
            for line in lines:
                chunks = line.split(" = ")
                #print("chunks",chunks)
                if (len(chunks) > 1):
                    chunk = chunks[1].strip()
                else:
                    chunk = None
                if ("mesh_style" in line):
                    self.config["mesh_style"] = chunk
                elif ("n_tacs_procs" in line):
                    self.config["n_tacs_procs"] = int(chunk)
                elif ("n_steps" in line):
                    self.config["nsteps"] = int(chunk)
                elif ("scenario_name" in line):
                    self.config["scenario_name"] = chunk
                elif ("csm_file" in line):
                    self.config["csmFile"] = chunk + ".csm"
                elif ("caps_constraint" in line):
                    self.config["capsConstraint"] = chunk
                elif ("constraint_type" in line):
                    self.config["constraintType"] = int(chunk)
                elif ("f2f_analysis" in line):
                    self.config["f2f_analysis_type"] = chunk
                elif ("fun3d_analysis" in line):
                    self.config["fun3d_analysis_type"] = chunk
                elif ("mach" in line):
                    self.config["mach"] = float(chunk)
                elif ("AOA" in line):
                    self.config["AOA"] = float(chunk)
                elif ("temperature" in line):
                    self.config["temperature"] = float(chunk)
                elif ("Re" in line):
                    self.config["Re"] = float(chunk)
                elif ("qinf" in line):
                    self.config["qinf"] = float(chunk)
                elif ("thermal_scale" in line):
                    self.config["thermal_scale"] = float(chunk)
                
            handle.close()


        #get the caps2fun environment directory
        self.caps2fun_dir = os.environ["CAPS2FUN"]
        self.src_dir = os.path.join(self.caps2fun_dir, "caps2fun")

        #read config on root proc
        if (self.comm.Get_rank() == 0):
            #initialize config attribute
            self.config = {}

            #read default config file
            default_cfg = os.path.join(self.src_dir, "default.cfg")
            parseFile(default_cfg)

            funtofem_folder = os.path.join(self.root_dir, "funtofem")
            if (not(os.path.exists(funtofem_folder))): os.mkdir(funtofem_folder)
            funtofem_cfg = os.path.join(funtofem_folder, "funtofem.cfg")
            if (not(os.path.exists(funtofem_cfg))): 
                shutil.copy(default_cfg, funtofem_cfg)
                sys.exit("The funtofem/funtofem.cfg file was missing. A copy has been created - edit and then run it again.\n")
            parseFile(funtofem_cfg)

        #now broadcast results from reading the config file
        self.config = self.comm.bcast(self.config, root=0)        

    def makeStatusFile(self):
        #initialize time
        self.start_time = time.time()

        #status file
        if (self.comm.Get_rank() == 0):
            statusFile = os.path.join(self.root_dir, "funtofem", "status.txt")
            self.status =  open(statusFile, "w")
            
        #into status
        self.cwrite("Running Funtofem with ESP/CAPS\n")

    def cwrite(self, text):
        if (self.comm.Get_rank() == 0):
            #write to the status file
            self.status.write(text)
        
            #immediately update it to be visibile in the file
            self.status.flush()

    def writeTime(self):
        dt = time.time() - self.start_time
        dt = round(dt)
        self.cwrite(", {} sec\n".format(dt))

    def commBarrier(self, location):
        print("proc #{} reached location {}\n".format(self.comm.Get_rank(), location))
        sys.stdout.flush()
        self.comm.Barrier()

    def readInput(self):
        #initialize DVdict as none for all procs
        self.DVdict = None
        self.functionNames = None
        self.mode = None
        self.eps = None
        self.x_dir = None

        #assume complex_step not used first
        self.complex = False

        if (self.comm.Get_rank() == 0):
            
            #print(curDir)
            funtofemFolder = os.path.join(self.root_dir, "funtofem")

            inputFile = os.path.join(funtofemFolder, "funtofem.in")
            inputHandle =  open(inputFile, "r")
            lines = inputHandle.readlines()

            self.functionNames = []

            #read in the DVdict from funtofem.in
            self.DVdict = []
            for line in lines:

                #read line if not empty
                if (len(line) > 1):
                    if ("function" in line):
                        parts = line.split(",")
                        name = parts[1]
                        name = name.strip()
                        self.functionNames.append(name)
                        ishape = 0
                        istruct = 0
                    elif ("mode" in line):
                        chunks = line.split(",")
                        mode = chunks[1].strip()
                        self.mode = mode
                        #if complex step mode, then turn on complex mode
                        if (self.mode == "complex_step"): 
                            self.complex = True
                            self.eps = 0
                            self.x_dir = []
                        else:
                            self.complex = False
                    elif ("eps" in line):

                        #read the epsilon from complex step run
                        chunks = line.split(",")
                        self.eps = float(chunks[1].strip())

                    elif ("x_dir" in line):

                        #read the x_direction for complex step run
                        chunks = line.split(",")
                        self.x_dir = np.zeros(len(chunks))
                        ct = 0
                        for chunk in chunks:
                            gatekeeper_chunk = "x_dir" in chunk
                            if (not(gatekeeper_chunk)):
                                self.x_dir[ct] = float(chunk)
                                ct += 1
                    else:
                        parts = line.split(",")
                        name = parts[0]
                        dvType = parts[1]
                        if (dvType == "shape"):
                            value = float(parts[2])
                            capsGroup = ""
                            DVind = ishape
                        elif (dvType == "struct"):
                            capsGroup = parts[2]
                            value = float(parts[3])
                            DVind = istruct
                        
                        #store DVdict
                        tempDict = {"name" : name,
                                    "type" : dvType,
                                    "capsGroup" : capsGroup,
                                    "value" : value,
                                    "ind" : DVind}
                        if (dvType == "shape"):
                            ishape += 1
                        elif (dvType == "struct"):
                            istruct += 1
                        self.DVdict.append(tempDict)


            inputHandle.close()

        #MPI broadcast from root proc
        self.DVdict = self.comm.bcast(self.DVdict, root=0)
        self.functionNames = self.comm.bcast(self.functionNames, root=0)
        self.mode = self.comm.bcast(self.mode, root=0)
        self.complex = self.comm.bcast(self.complex,root=0)

        if (self.mode == "adjoint"):
            #self.cwrite("Running in real mode\n")
            pass

        if (self.mode == "complex_step"):
            #self.cwrite("Running in complex mode\n")

            self.eps = self.comm.bcast(self.eps, root=0)
            self.x_dir = self.comm.bcast(self.x_dir,root=0)

    def initializeAIMs(self):
        if (self.comm.Get_rank() == 0):
            #initialize all 6 ESP/CAPS AIMs used for the fluid and structural analysis

            #initialize pyCAPS structural problem
            self.capsStruct = pyCAPS.Problem(problemName = "CAPS_struct",
                        capsFile = self.config["csmFile"],
                        outLevel = 1)

            self.capsStruct.geometry.cfgpmtr["cfdOn"].value = 0

            #print the structure mesh deskeys
            print("Design keys... {}".format(self.capsStruct.geometry.despmtr.keys()))

            #self.cwrite("Initialized caps Struct AIM\n")

            #initialize pyCAPS fluid problem
            self.capsFluid = pyCAPS.Problem(problemName = "CAPS_fluid",
                        capsFile = self.config["csmFile"],
                        outLevel = 1)
            self.capsFluid.geometry.cfgpmtr["cfdOn"].value = 1
            #self.cwrite("Initialized caps fluid AIM\n")

            #initialize egads Aim
            self.egadsAim = self.capsStruct.analysis.create(aim="egadsTessAIM")
            #self.cwrite("Initialized egads AIM\n")

            #initialize tacs Aim
            self.tacsAim = self.capsStruct.analysis.create(aim = "tacsAIM", name = "tacs")
            #self.cwrite("Initialized tacs AIM\n")


            if (self.config["mesh_style"] == "pointwise"):
                #initialize pointwise AIM
                self.pointwiseAim = self.capsFluid.analysis.create(aim = "pointwiseAIM",
                                            name = "pointwise")
                #self.cwrite("Initialized pointwise AIM\n")
            elif (self.config["mesh_style"] == "tetgen"):
                #initialize egads fluid aim
                self.egadsFluidAim = self.capsFluid.analysis.create(aim = "egadsTessAIM", name = "egadsTess")
                #self.cwrite("Initialized egadsTessAIM for fluid\n")
                
                #initialize tetgen aim
                self.tetgenAim = self.capsFluid.analysis.create(aim = "tetgenAIM", name = "tetgen")
                #self.cwrite("Initialized tetgen AIM\n")

            #initialize FUN3D AIM from Pointwise mesh
            self.fun3dAim = self.capsFluid.analysis.create(aim = "fun3dAIM",
                                    name = "fun3d")
            #self.cwrite("Initialized fun3d AIM\n")

            #structural mesh settings
            self.structureMeshSettings()
            #self.cwrite("Set Structure mesh settings\n")

            #fluid mesh settings
            self.fluidMeshSettings()
            #self.cwrite("Set Fluid mesh settings\n")

            #set fun3d settings
            self.fun3dSettings()
            #self.cwrite("Set fun3d settings\n")


    def structureMeshSettings(self):
        #names the bdf and dat files as pointwise.ext or tetgen.ext
        self.tacsAim.input.Proj_Name = self.config["mesh_style"]
        
        #Egads Aim section, for mesh
        self.egadsAim.input.Edge_Point_Min = 5
        self.egadsAim.input.Edge_Point_Max = 20

        self.egadsAim.input.Mesh_Elements = "Quad"

        self.egadsAim.input.Tess_Params = [0.25,.01,15]

        #increase the precision in the BDF file
        self.tacsAim.input.File_Format = "Large"
        self.tacsAim.input.Mesh_File_Format = "Large"

        # Link the mesh
        self.tacsAim.input["Mesh"].link(self.egadsAim.output["Surface_Mesh"])

        # Set analysis type
        self.tacsAim.input.Analysis_Type = "Static"

        #materials section    
        aluminum    = {"materialType" : "isotropic",
                        "youngModulus" : 72.0E9 ,
                        "poissonRatio": 0.33,
                        "density" : 2.8E3,
                        "tensionAllow" :  20.0e7}

        self.tacsAim.input.Material = {"aluminum": aluminum}

        # Material properties section
        propDict = {}
        for DV in self.DVdict:
            if (DV["type"] == "struct"):
                capsGroup = DV["capsGroup"]
                if (len(capsGroup) > 0):

                    bendingInertiaRatio = 1.0 #default
                    shearMembraneRatio = 5.0/6.0 #default

                    #factor to artificially boost OML bending and membrane stiffnesses as if stringers were there
                    boostFactor = 5.0

                    #artificial boost to bending inertia as if stringers were there in OML
                    if ("OML" in capsGroup): 
                        bendingInertiaRatio *= boostFactor
                        shearMembraneRatio *= boostFactor

                        print("applied boost factor to capsGroup {}".format(capsGroup))

                    #make the shell property
                    shell = {"propertyType" : "Shell",
                        "membraneThickness" : DV["value"],
                        "material"        : "aluminum",
                        "bendingInertiaRatio" : bendingInertiaRatio, # Default
                        "shearMembraneRatio"  : shearMembraneRatio} # Default

                    
                    propDict[capsGroup] = shell

        self.tacsAim.input.Property = propDict

        # constraint section
        constraint1 = {"groupName" : self.config["capsConstraint"],
                        "dofConstraint" : self.config["constraintType"]} #123

        self.tacsAim.input.Constraint = {"fixRoot": constraint1}

        # don't set any loads

        ##-----------setup DVs and DVRs-----------##

        #make initial DV and DVR dict
        DVdict = {}
        DVRdict = {}

        #add thickDVs and shapeDVs to caps
        for DV in self.DVdict:
            dvname = DV["name"]
            if (DV["type"] == "struct"):
                #add thickDV entry into DV_Relations and DV Dicts
                DVRdict[dvname] = self.makeThicknessDVR(dvname)
                DVdict[dvname] = self.makeThicknessDV(DV["capsGroup"], DV["value"])
            elif (DV["type"] == "shape"): #geomDV, add empty entry into DV dicts
                DVdict[dvname] = {}

        #input DVdict and DVRdict into tacsAim
        self.tacsAim.input.Design_Variable = DVdict
        self.tacsAim.input.Design_Variable_Relation = DVRdict

    def fluidMeshSettings(self):
        if (self.config["mesh_style"] == "pointwise"):

            #wall bc settings (wall is the OML)
            if (self.config["fun3d_analysis_type"] == "inviscid"):
                self.wallBC = {"bcType" : "inviscid"}
                wallSpacing = 0.1
            elif (self.config["fun3d_analysis_type"] == "laminar"):
                wallSpacing = 0.01
                self.wallBC = {"bcType" : "viscous",
                "boundaryLayerSpacing" : wallSpacing}
            elif (self.config["fun3d_analysis_type"] == "turbulent"):
                wallSpacing = 0.01
                self.wallBC = {"bcType" : "viscous",
                "boundaryLayerSpacing" : wallSpacing}

            # Dump VTK files for visualization
            self.pointwiseAim.input.Proj_Name   = "TransportWing"
            self.pointwiseAim.input.Mesh_Format = "VTK"

            # Connector level
            self.pointwiseAim.input.Connector_Turn_Angle       = 1
            self.pointwiseAim.input.Connector_Prox_Growth_Rate = 1.2
            self.pointwiseAim.input.Connector_Source_Spacing   = True

            # Domain level
            self.pointwiseAim.input.Domain_Algorithm    = "AdvancingFront"
            self.pointwiseAim.input.Domain_Max_Layers   = 15
            self.pointwiseAim.input.Domain_Growth_Rate  = 1.75
            self.pointwiseAim.input.Domain_TRex_ARLimit = 40.0 #def 40.0, lower inc mesh size
            self.pointwiseAim.input.Domain_Decay        = 0.5
            self.pointwiseAim.input.Domain_Iso_Type = "Triangle" #"TriangleQuad"
            self.pointwiseAim.input.Domain_Wall_Spacing = wallSpacing

            # Block level
            self.pointwiseAim.input.Block_Boundary_Decay       = 0.5
            self.pointwiseAim.input.Block_Collision_Buffer     = 1.0
            self.pointwiseAim.input.Block_Max_Skew_Angle       = 170.0
            self.pointwiseAim.input.Block_Edge_Max_Growth_Rate = 2.0
            self.pointwiseAim.input.Block_Full_Layers          = 1
            self.pointwiseAim.input.Block_Max_Layers           = 100
            self.pointwiseAim.input.Block_TRexType = "TetPyramid"
            #T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms)        

        
            self.pointwiseAim.input.Mesh_Sizing = {"wall": self.wallBC,
                "Farfield": {"bcType":"Farfield"}}

        elif (self.config["mesh_style"] == "tetgen"):
            self.egadsFluidAim.input.Tess_Params = [0.01, 0.01, 0]
            self.tetgenAim.input.Quality_Rad_Edge = 1.0
            self.egadsFluidAim.input.Mesh_Sizing = {"Farfield": {"tessParams" : [0.005, 0.01, 0]}}
            self.egadsFluidAim.input.Mesh_Sizing = {"Symmetry": {"tessParams" : [0.0005, 0.01, 0]}}
            self.egadsFluidAim.input.Mesh_Sizing = {"wall": {"tessParams" : [0.001, 0.01, 0]}}

            #tetgen solvers['flow'] = fun3d_interface(self.comm, )AIM mesh params
            self.tetgenAim.input.Preserve_Surf_Mesh = True
            self.tetgenAim.input["Surface_Mesh"].link(self.egadsFluidAim.output["Surface_Mesh"])
            self.tetgenAim.input.Mesh_Format = "AFLR3"

    def fun3dSettings(self):

        
        self.fun3dAim.input.Boundary_Condition = {"wall": self.wallBC,
                "Farfield": {"bcType":"Farfield"}}

        #add thickDVs and geomDVs to caps
        DVdict = {}
        for DV in self.DVdict:
            dvname = DV["name"]
            if (DV["type"] == "shape"): #geomDV, add empty entry into DV dicts
                DVdict[dvname] = {}

        #input DVdict and DVRdict into tacsAim
        if (len(DVdict) > 0): self.fun3dAim.input.Design_Variable = DVdict

        #fun3d design sensitivities settings
        self.fun3dAim.input.Design_SensFile = True
        self.fun3dAim.input.Design_Sensitivity = True

        #############################
        # namelist and mapbc settings
        self.capsFluid.analysis["fun3d"].input.Overwrite_NML = False
        self.fun3dnml = f90nml.Namelist()
        
        #project section
        self.fun3dnml["project"] = f90nml.Namelist()
        self.fun3dnml["project"]["project_rootname"] = self.fun3dAim.input.Proj_Name
        #self.fun3dAim.input.Proj_Name = self.config["mesh_style"]

        #governing equation section
        self.fun3dnml["governing_equations"] = f90nml.Namelist()

        #apply fun3d analysis type
        #options = "inviscid", "laminar", "turbulent"
        if (self.config["fun3d_analysis_type"] == "inviscid"):
            self.fun3dnml["governing_equations"]["viscous_terms"] = "inviscid"
        elif (self.config["fun3d_analysis_type"] == "laminar"):
            self.fun3dnml["governing_equations"]["eqn_type"] = "compressible"
            self.fun3dnml["governing_equations"]["viscous_terms"] = "laminar"
        elif (self.config["fun3d_analysis_type"] == "turbulent"):
            pass

        #raw grid section
        self.fun3dnml["raw_grid"] = f90nml.Namelist()
        #self.fun3dAim.input.Mesh_ASCII_Flag = False
        self.fun3dnml["raw_grid"]["grid_format"] = "aflr3"
        self.fun3dnml["raw_grid"]["data_format"] = "default"
        self.fun3dnml["raw_grid"]["swap_yz_axes"] = False

        #reference physical properties section
        self.fun3dnml["reference_physical_properties"] = f90nml.Namelist()
        self.fun3dnml["reference_physical_properties"]["mach_number"] = self.config["mach"]
        self.fun3dnml["reference_physical_properties"]["angle_of_attack"] = self.config["AOA"]
        self.fun3dnml["reference_physical_properties"]["reynolds_number"] = self.config["Re"]
        self.fun3dnml["reference_physical_properties"]["temperature"] = self.config["temperature"]
        self.fun3dnml["reference_physical_properties"]["temperature_units"] = "Kelvin"
        #self.capsFluid.analysis["fun3d"].input.Alpha = 1.0
        #self.capsFluid.analysis["fun3d"].input.Mach = 0.5
        #self.capsFluid.analysis["fun3d"].input.Re = 35e6
        
        #inviscid flux method section
        self.fun3dnml["inviscid_flux_method"] = f90nml.Namelist()
        self.fun3dnml["inviscid_flux_method"]["flux_construction"] = "roe"
        self.fun3dnml["inviscid_flux_method"]["flux_limiter"] = "hminmod"
        self.fun3dnml["inviscid_flux_method"]["smooth_limiter_coeff"] = 1.0
        self.fun3dnml["inviscid_flux_method"]["freeze_limiter_iteration"] = int(5.0/6 * self.config["nsteps"])

        #nonlinear solver parameters section
        self.fun3dnml["nonlinear_solver_parameters"] = f90nml.Namelist()
        self.fun3dnml["nonlinear_solver_parameters"]["schedule_iteration"] = [1, 80]
        self.fun3dnml["nonlinear_solver_parameters"]["schedule_cfl"] = [2, 100]
        if (not(self.config["fun3d_analysis_type"] == "inviscid")):
            self.fun3dnml["nonlinear_solver_parameters"]["time_accuracy"] = "steady"
            self.fun3dnml["nonlinear_solver_parameters"]["time_step_nondim"] = 0.1
            self.fun3dnml["nonlinear_solver_parameters"]["subiterations"] = 0
            self.fun3dnml["nonlinear_solver_parameters"]["schedule_cflturb"] = [50.0,50.0]


        #force moment integ properties section
        self.fun3dnml["force_moment_integ_properties"] = f90nml.Namelist()
        self.fun3dnml["force_moment_integ_properties"]["area_reference"] = self.capsStruct.geometry.despmtr["area"].value / 2

        #code run and control section
        self.fun3dnml["code_run_control"] = f90nml.Namelist()
        self.fun3dnml["code_run_control"]["steps"] = self.config["nsteps"]
        self.fun3dnml["code_run_control"]["stopping_tolerance"] = 1.0e-15
        self.fun3dnml["code_run_control"]["restart_write_freq"] = 1000
        self.fun3dnml["code_run_control"]["restart_read"] = "off"

        #global settings
        self.fun3dnml["global"] = f90nml.Namelist()
        self.fun3dnml["global"]["moving_grid"] = True
        self.fun3dnml["global"]["volume_animation_freq"] = -1
        self.fun3dnml["global"]["boundary_animation_freq"] = -1

        #mesh elasticity settings
        self.fun3dnml["elasticity_gmres"] = f90nml.Namelist()
        self.fun3dnml["elasticity_gmres"]["algebraic_mesh_deform"] = False #default False
        self.fun3dnml["elasticity_gmres"]["nsearch"] = 200 #default 50
        self.fun3dnml["elasticity_gmres"]["tol"] = 1.e-14
        self.fun3dnml["elasticity_gmres"]["deformation_substeps"] = 1 #default 1, other value 5
        self.fun3dnml["elasticity_gmres"]["deform_from_initial_mesh"] = True #default true
        self.fun3dnml["elasticity_gmres"]["use_substeps_each_step"] = False #default False, need to turn to True if deformation_substeps used
        self.fun3dnml["elasticity_gmres"]["elasticity"] = 1 #default 1, option 2
        self.fun3dnml["elasticity_gmres"]["elasticity_exponent"] = 1.0 #default 1.0, change to 2.0 if needed
        self.fun3dnml["elasticity_gmres"]["nrestarts"] = 1
        self.fun3dnml["elasticity_gmres"]["poisson_ratio"] = 0.0 #default 0.0
        
        #massoud output settings
        self.fun3dnml["massoud_output"] = f90nml.Namelist()
        self.fun3dnml["massoud_output"]["funtofem_include_skin_friction"] = False

        #volume output variables
        self.fun3dnml["volume_output_variables"] = f90nml.Namelist()
        self.fun3dnml["volume_output_variables"]["export_to"] = "vtk"
        self.fun3dnml["volume_output_variables"]["x"] = False
        self.fun3dnml["volume_output_variables"]["y"] = False
        self.fun3dnml["volume_output_variables"]["z"] = False
        self.fun3dnml["volume_output_variables"]["temperature"] = True
        self.fun3dnml["volume_output_variables"]["mach"] = True
        self.fun3dnml["volume_output_variables"]["p"] = True

        #boundary output variables
        self.fun3dnml["boundary_output_variables"] = f90nml.Namelist()
        #boundary list indexes probably auto set from fun3dAim
        self.fun3dnml["boundary_output_variables"]["number_of_boundaries"] = -1
        self.fun3dnml["boundary_output_variables"]["boundary_list"] = "1-2"
        self.fun3dnml["boundary_output_variables"]["temperature"] = True
        self.fun3dnml["boundary_output_variables"]["mach"] = True
        self.fun3dnml["boundary_output_variables"]["p"] = True

        ##############################
        # fun3d settings for moving_body.input file
        self.moving_body_input = f90nml.Namelist()

        #moving body settings for funtofem to fun3d
        bodyName = self.config["csmFile"].split(".")[0]
        nBodies = 1
        nBoundaries = 1
        bndryArray = [[2]]
        bndryArray = list(bndryArray)

        #body definitions
        self.moving_body_input["body_definitions"] = f90nml.Namelist()
        self.moving_body_input["body_definitions"]["n_moving_bodies"] = nBodies
        self.moving_body_input["body_definitions"]["body_name"] = [bodyName]
        self.moving_body_input["body_definitions"]["parent_name"] = [""] # '' means motion relative to inertial ref frame
        self.moving_body_input["body_definitions"]["n_defining_bndry"] = [nBoundaries] #number of boundaries that define this body
        self.moving_body_input["body_definitions"]["defining_bndry(1,1)"] = 2 #index 1: boundary number index 2: body number
        self.moving_body_input["body_definitions"]["motion_driver"] = ["funtofem"] #tells fun3d to use motion inputs from python
        self.moving_body_input["body_definitions"]["mesh_movement"] = ["deform"] #can use 'rigid', 'deform', 'rigid+deform' with funtofem interface

        ##############################
        # fun3d settings for moving_body.input file
        self.moving_body_input = f90nml.Namelist()

        #moving body settings for funtofem to fun3d
        bodyName = self.config["csmFile"].split(".")[0]
        nBodies = 1
        nBoundaries = 1
        bndryArray = [[2]]
        bndryArray = list(bndryArray)

        #body definitions
        self.moving_body_input["body_definitions"] = f90nml.Namelist()
        self.moving_body_input["body_definitions"]["n_moving_bodies"] = nBodies
        self.moving_body_input["body_definitions"]["body_name"] = [bodyName]
        self.moving_body_input["body_definitions"]["parent_name"] = [""] # '' means motion relative to inertial ref frame
        self.moving_body_input["body_definitions"]["n_defining_bndry"] = [nBoundaries] #number of boundaries that define this body
        self.moving_body_input["body_definitions"]["defining_bndry(1,1)"] = 2 #index 1: boundary number index 2: body number
        self.moving_body_input["body_definitions"]["motion_driver"] = ["funtofem"] #tells fun3d to use motion inputs from python
        self.moving_body_input["body_definitions"]["mesh_movement"] = ["deform"] #can use 'rigid', 'deform', 'rigid+deform' with funtofem interface

    def initF2F(self):
        #sort struct DVs to match alphabetic, numeric sorting of ESP/CAPS
        structDVs = self.sortStructDVs()
        
        maximum_mass = 40.0 
        num_tacs_dvs = len(structDVs)

        # Set up the communicators
        n_procs = self.comm.Get_size()
        if (self.config["n_tacs_procs"] > n_procs): self.config["n_tacs_procs"] = n_procs

        world_rank = self.comm.Get_rank()
        if world_rank < self.config["n_tacs_procs"]:
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
        self.wing = Body('wing', analysis_type=self.config["f2f_analysis_type"], group=0,boundary=2)

        for i in range(num_tacs_dvs):
            self.wing.add_variable('structural',Variable('thickness '+ str(i),value=structDVs[i],lower = 0.0001, upper = 1.0))

        self.model.add_body(self.wing)

        myscenario = Scenario(self.config["scenario_name"], group=0, steady=True, steps=self.config["nsteps"])

        #add functions to scenario
        for functionName in self.functionNames:
            #strip newline character from function name
            functionName = functionName.strip()

            if (functionName in ["cl","cd"]):
                function = Function(functionName,analysis_type='aerodynamic')
            elif (functionName in ["ksfailure"]):
                function = Function(functionName,analysis_type='structural')
            elif (functionName in ["mass"]):
                function = Function(functionName,analysis_type='structural',adjoint=False)
            else:
                sys.exit("Function name {} not in accepted functions".format(functionName))
            myscenario.add_function(function)

        self.model.add_scenario(myscenario)

        #==================================================================================================#

        #broadcast directories and datFile to rest of procs
        fun3d_parent_dir = None
        datFile = None
        if (self.comm.Get_rank() == 0):
            fun3d_parent_dir = os.path.join(self.fun3dAim.analysisDir, "..")
            datFile = os.path.join(self.tacsAim.analysisDir, self.config["mesh_style"] + ".dat")
        fun3d_parent_dir = self.comm.bcast(fun3d_parent_dir, root=0)
        datFile = self.comm.bcast(datFile, root=0)

        solvers = {}
        solvers['flow'] = Fun3dInterface(self.comm,self.model,flow_dt=1.0, qinf=self.config["qinf"], thermal_scale=self.config["thermal_scale"], fun3d_dir=fun3d_parent_dir)
        solvers['structural'] = TACSinterface(self.comm,tacs_comm,self.model,self.config["n_tacs_procs"], datFile, structDVs)

        # L&D transfer options
        transfer_options = {'analysis_type': self.config["f2f_analysis_type"],
                            'scheme': 'meld', 'thermal_scheme': 'meld'}

        # instantiate the driver
        self.driver = FUNtoFEMnlbgs(solvers,self.comm,tacs_comm,0,self.comm,0,transfer_options,model=self.model)

        #section to update funtofem DV values for next fun3d run
        self.model.set_variables(structDVs)

    def sortStructDVs(self):
        def capsCompare(elem):
            return elem["capsGroup"]


        structDVnames = []
        #make list of DV values and names
        for DV in self.DVdict:
            if (DV["type"] == "struct"):
                structDVnames.append(DV["name"])
        
        #sort the names and values based on ESP/CAPS sorting
        self.DVdict.sort(key=capsCompare)

        #get the sorted structDVs
        structDVs = []
        inds = []
        for DV in self.DVdict:
            if (DV["type"] == "struct"):
                structDVs.append(DV["value"])
                inds.append(DV["ind"])

        #if complex step mode add perturbation
        if (self.complex):
            for i in range(len(structDVs)):
                #ith struct DV has original index ind, thus from x_dir
                structDVs[i] += 1j * self.eps * self.x_dir[inds[i]]

        return structDVs

    def forwardAnalysis(self):
        #update status
        self.cwrite("----------------------------")
        self.cwrite("----------------------------\n")

        #set design variables from desvarDict
        self.updateDesign()

        #generate structure mesh with egads and tacs AIMs
        self.start_time = time.time()
        self.buildStructureMesh()
        self.cwrite("built structure mesh")
        self.writeTime()

        #generate fluid mesh with Pointwise
        self.start_time = time.time()

        self.buildFluidMesh()
        self.cwrite("built fluid mesh")
        self.writeTime()


        self.start_time = time.time()
        #prepare fun3d config files such as mapbc, nml
        self.runFun3dConfig()        

        #initialize fun2fem with new meshes
        self.initF2F()

        #run FUNtoFEM forward analysis
        self.cwrite("Running F2F forward analysis... ")
        self.driver.solve_forward()
        self.functions = self.model.get_functions()

        self.cwrite("completed F2F forward analysis")
        dt = time.time() - self.start_time
        dtPerStep = round(dt/self.config["nsteps"])
        self.cwrite(", {} sec/step".format(dtPerStep))
        self.writeTime()

    def adjointAnalysis(self):
        #update status
        self.start_time = time.time()
        self.cwrite("Running F2F adjoint analysis... ")

        #run funtofem adjoint analysis, with fun3d
        self.driver.solve_adjoint()

        #update status
        self.cwrite("finished adjoint analysis")
        dt = time.time() - self.start_time
        dtPerStep = round(dt/self.config["nsteps"])
        self.cwrite(", {} sec/step".format(dtPerStep))
        self.writeTime()

        #get the funtofem gradients
        f2fgrads = self.model.get_function_gradients()
        self.start_time = time.time()

        #number of shape DVs
        self.nshapeDV = 0
        self.nstructDV = 0
        for DV in self.DVdict:
            if (DV["type"] == "shape"): self.nshapeDV += 1
            if (DV["type"] == "struct"): self.nstructDV += 1

        #update gradient for non shape DVs
        nfunc = len(self.functions)
        self.structGrad = np.zeros((nfunc, self.nstructDV))
        
        for funcInd in range(nfunc):
            ct = 0
            for DV in self.DVdict:
                dvname = DV["name"]
                dvind = DV["ind"]
                if (not(DV["type"] == "shape")): 
                    self.structGrad[funcInd, dvind] = f2fgrads[funcInd][ct].real
                    ct += 1        

        if (self.nshapeDV > 0):
            #get aero and struct mesh sensitivities for shape DVs
            self.aeroIds, self.aero_mesh_sens = self.wing.collect_coordinate_derivatives(self.comm, "aero")
            self.structIds, self.struct_mesh_sens = self.wing.collect_coordinate_derivatives(self.comm, "struct")

            self.commBarrier("collected coordinate derivatives")

            if (self.comm.Get_rank() == 0):

                #compute shape derivatives from aero and struct mesh sensitivities
                self.shapeGrad = self.computeShapeDerivatives()
                self.cwrite("computed shape derivative chain rule products")
                self.writeTime()

            else:
                self.shapeGrad = None

            #barrier to wait for root proc to finish computing the shape derivatives
            self.commBarrier("computeShapeDerivatives()")

            self.shapeGrad = self.comm.bcast(self.shapeGrad, root=0)
        
        else: #no shape DVs
            self.nfunc = len(self.functions)
            self.shapeGrad = np.zeros((self.nfunc, 0))

        

        #send or store gradient
        #return structGrad, self.shapeGrad

    def updateDesign(self):
        
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

                    #for aim in [self.capsStruct, self.capsFluid, self.pointwiseAim, self.fun3dAim, self.tacsAim]:
                    for aim in [self.capsStruct, self.capsFluid]:
                        aim.geometry.despmtr[dvname].value = value

                elif (DV["type"] == "struct"):
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

    def buildFluidMesh(self):
        self.cwrite("Building fluid mesh... ")

        if (self.comm.Get_rank() == 0):
            if (self.config["mesh_style"] == "pointwise"):
                #build fluid mesh by running pointwise and then linking with fun3d
                self.runPointwise()
                self.cwrite("ran pointwise, ")

            #update fun3d with current mesh so it knows mesh sensitivity
            if (self.config["mesh_style"] == "pointwise"):
                self.fun3dAim.input["Mesh"].link(self.pointwiseAim.output["Volume_Mesh"])
            elif (self.config["mesh_style"] == "tetgen"):
                self.fun3dAim.input["Mesh"].link(self.tetgenAim.output["Volume_Mesh"])
            self.fun3dAim.preAnalysis()
            self.cwrite("linked to fun3d, ")

    def runFun3dConfig(self):
        #build fun3d config files, mapbc and nml
        if (self.comm.Get_rank() == 0):
            #set caps and funtofem flow folders for fun3d
            caps_flow_dir = os.path.join(self.fun3dAim.analysisDir, "Flow")

            #run the namelist generator, and fun3d aim preanalysis to build fun3d.nml file
            self.fun3dnml.write(os.path.join(caps_flow_dir, "fun3d.nml"), force=True)

            #write the moving_body.input file
            self.moving_body_input.write(os.path.join(caps_flow_dir, "moving_body.input"), force=True)

            #add names of the BCs to the mapbc file
            mapbc = os.path.join(caps_flow_dir, "fun3d_CAPS.mapbc")
            mapbc_hdl = open(mapbc, "r")
            lines = mapbc_hdl.readlines()
            mapbc_hdl.close()
            wr_hdl = open(mapbc,"w")
            newlines = []
            bc_inds = ["1","2"]
            bc_names = ["Farfield", "wall"]
            ind = 0
            for line in lines:
                chunks = line.split(" ")
                if (len(chunks) > 1):
                    line = line.strip()
                    if (bc_inds[ind] in chunks[0]):
                        line += " " + bc_names[ind]
                        ind += 1
                    line += "\n"

                wr_hdl.write(line)
            wr_hdl.close()
                    
            #get the caps2fun project dir
            caps2fun_proj_dir = os.environ["CAPS2FUN"]
            archive_folder = os.path.join(caps2fun_proj_dir, "archive")

            #move perturb.input from archive folder if complex mode
            if (self.complex):
                src = os.path.join(archive_folder,"perturb.input")
                dest = os.path.join(caps_flow_dir, "perturb.input")
                shutil.copy(src, dest)


    def runPointwise(self):
        #run AIM pre-analysis
        self.pointwiseAim.preAnalysis()

        #move to test directory
        os.chdir(self.pointwiseAim.analysisDir)

        CAPS_GLYPH = os.environ["CAPS_GLYPH"]
        for i in range(1): #can run extra times if having license issues
                os.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")
            #if os.path.isfile('caps.GeomToMesh.gma') and os.path.isfile('caps.GeomToMesh.ugrid'): break

        os.chdir(self.root_dir)

        #run AIM postanalysis, files in self.pointwiseAim.analysisDir
        self.pointwiseAim.postAnalysis()      


    def runF2F(self):

        #run forward analysis
        self.forwardAnalysis()        

         #count the number of design variables and functions
        self.countDV()

        #run adjoint analysis
        if (self.mode == "adjoint"):
            self.adjointAnalysis()

        #write to output
        self.writeOutput()

        #write that you finished the run
        self.cwrite("Finished the funtofem call\n")

        #close status file
        if (self.comm.Get_rank() == 0): self.status.close()

    def writeOutput(self):
        #write output file, funtofem.out
        #for modes: adjoint, complex_step, or forward

        if (self.comm.Get_rank() == 0):

            #...write the output functions, gradients, etc
            funtofemFolder = os.path.join(self.root_dir, "funtofem")

            outputFile = os.path.join(funtofemFolder, "funtofem.out")

            outputHandle =  open(outputFile, "w")

            if (self.mode == "adjoint"):

                #format of the output file:
                #nfunc,nDV
                #func,name,value
                #grad,DVname,deriv1
                #grad,DVname,deriv2
                #grad,DVname,deriv3
                #...
                #grad,DVname,derivN
                #... for each function

                #write number of functions, variables, etc.
                nDV = self.nshapeDV + self.nstructDV
                outputHandle.write("{},{}\n".format(self.nfunc, nDV))
                ifunc = 0
                for func in self.functions:

                    name = self.functionNames[ifunc]
                    value = self.functions[ifunc].value.real

                    #write function name and value
                    outputHandle.write("func,{},{}\n".format(name,value))

                    #write shape gradient
                    for ishape in range(self.nshapeDV):
                        #find that shape variable
                        for DV in self.DVdict:
                            isShape = DV["type"] == "shape"
                            matchingInd = DV["ind"] == ishape
                            if (isShape and matchingInd):
                                name = DV["name"]
                                deriv = self.shapeGrad[ifunc, ishape]
                                outputHandle.write("grad,{},{}\n".format(name,deriv))
                    #write struct gradient
                    for istruct in range(self.nstructDV):
                        #find each struct variable if out of order
                        for DV in self.DVdict:
                            isStruct = DV["type"] == "struct"
                            matchingInd = DV["ind"] == istruct
                            if (isStruct and matchingInd):
                                name = DV["name"]
                                deriv = self.structGrad[ifunc, istruct]
                                outputHandle.write("grad,{},{}\n".format(name,deriv))

                    #update function counter
                    ifunc += 1

            elif (self.mode == "forward"):
                #format of the output file:
                #nfunc,nDV
                #func,name,value
                #write number of functions, variables, etc.
                nDV = self.nshapeDV + self.nstructDV
                outputHandle.write("{},{}\n".format(self.nfunc, nDV))
                ifunc = 0
                for func in self.functions:

                    name = self.functionNames[ifunc]
                    value = self.functions[ifunc].value.real

                    #write function name and value
                    outputHandle.write("func,{},{}\n".format(name,value))

                    #update function counter
                    ifunc += 1


            elif (self.mode == "complex_step"):

                #write the following
                #nfunc
                #func,fname,real,value,imag,value
                #for each function

                #write the number of functions
                nfunc = len(self.functions)
                outputHandle.write("{}\n".format(nfunc))

                #write the value of each function
                ifunc = 0
                for function in self.functions:

                    #get function name
                    name = self.functionNames[ifunc]
                    ifunc += 1

                    #get real and imaginary parts, with imaginary part normalized by epsilon
                    real_part = function.value.real
                    imag_part = function.value.imag/self.eps

                    #write the line to the output file
                    line = "func,{},real,{},imag,{}\n".format(name,real_part,imag_part)
                    outputHandle.write(line)

            #close the output file
            outputHandle.close()

    def computeShapeDerivatives(self):
        self.cwrite("Compute shape derivatives... ")

        #initialize shape gradient again at zero
        self.initShapeGrad()

        if (self.comm.Get_rank() == 0):

            #add struct_mesh_sens part to shape DV derivatives#struct shape derivatives
            self.applyStructMeshSens()

            #add aero_mesh_sens part to shape DV derivatives
            self.applyAeroMeshSens()

        return self.shapeGrad

    def countDV(self):
        self.nfunc = len(self.functions)

        #determine number of shapeDV
        self.nshapeDV = 0
        self.nstructDV = 0
        for DV in self.DVdict:
            if (DV["type"] == "shape"): self.nshapeDV += 1
            if (DV["type"] == "struct"): self.nstructDV += 1

    def initShapeGrad(self):
        self.countDV()

        self.shapeGrad = np.zeros((self.nfunc, self.nshapeDV))

        self.cwrite("initialized shape gradient\n")

    def applyStructMeshSens(self):
        #print struct mesh sens to struct mesh sens file
        #where to print .sens file
        structSensFile = os.path.join(self.tacsAim.analysisDir, self.tacsAim.input.Proj_Name+".sens")

        #open the file
        with open(structSensFile, "w") as f:
            
            #write (nfunctions) in first line
            f.write("{}\n".format(self.nfunc))
            
            #for each function mass, stress, etc.
            for funcInd in range(self.nfunc):
                
                #get the pytacs/tacs sensitivity w.r.t. mesh for that function
                sens = self.struct_mesh_sens

                #write the key,value,nnodes of the function
                funcKey = "func#" + str(funcInd)
                f.write(funcKey + "\n")
                f.write("{}\n".format(self.functions[funcInd].value.real))
                f.write("{}\n".format(len(self.structIds)))

                #for each node, print nodeind, dfdx, dfdy, dfdz for that mesh element
                ct = 0
                for nodeind in self.structIds: # d(Func1)/d(xyz)
                    bdfind = nodeind + 1
                    f.write("{} {} {} {}\n".format(bdfind, sens[3*ct, funcInd].real, sens[3*ct+1, funcInd].real, sens[3*ct+2, funcInd].real))
                    ct += 1
        
            f.close()

        #update status
        self.cwrite("\tprinted struct.sens file, ")

        #run aim postanalysis
        self.tacsAim.postAnalysis()
        self.cwrite("completed tacsAim postAnalysis(), ")

        #update shape DV derivatives from struct mesh part
        for funcInd in range(self.nfunc):
            funcKey = "func#" + str(funcInd)
            dvct = 0
            for DV in self.DVdict:
                dvname = DV["name"]
                if (DV["type"] == "shape"): 
                    self.shapeGrad[funcInd, dvct] += self.tacsAim.dynout[funcKey].deriv(dvname)
                    dvct += 1

        #update status
        self.cwrite("finished struct mesh contribution to shape DVs\n")


    def applyAeroMeshSens(self):
        #print aero mesh sens to aero mesh sens file
        #where to print .sens file
        aeroSensFile = os.path.join(self.fun3dAim.analysisDir, self.fun3dAim.input.Proj_Name+".sens")
        #print("Writing aero sens file, {}".format(aeroSensFile))

        #make aero mesh derivatives (surface aero mesh) one based
        aero_nnodes = len(self.aeroIds)
        #minId = min(self.aeroIds)
        #maxId = max(self.aeroIds)
        
        #open the file
        with open(aeroSensFile, "w") as f:
            
            #write (nfunctions) in first line
            f.write("{}\n".format(self.nfunc))
            
            #for each function mass, stress, etc.
            for funcInd in range(self.nfunc):
                
                #get the pytacs/tacs sensitivity w.r.t. mesh for that function
                sens = self.aero_mesh_sens

                #write the key,value,nnodes of the function
                funcKey = "func#" + str(funcInd)
                f.write(funcKey + "\n")
                f.write("{}\n".format(self.functions[funcInd].value.real))
                
                #print aero_nnodes on aerodynamic surface mesh
                aero_nnodes = len(self.aeroIds)
                f.write("{}\n".format(aero_nnodes))

                #for each node, print nodeind, dfdx, dfdy, dfdz for that mesh element
                ct = 0
                for nodeind in self.aeroIds: # d(Func1)/d(xyz)
                    #bdfind = nodeind - (minId - 1)
                    bdfind = nodeind
                    f.write("{} {} {} {}\n".format(bdfind, sens[3*ct, funcInd].real, sens[3*ct+1, funcInd].real, sens[3*ct+2, funcInd].real))
                    ct += 1

            f.close()

        #update status
        self.cwrite("\tprinted aero.sens file, ")

        #run aim postanalysis
        self.fun3dAim.postAnalysis()
        self.cwrite("completed pointwiseAim postAnalysis(), ")

        #update shape DV derivatives from aero mesh part
        for funcInd in range(self.nfunc):
            funcKey = "func#" + str(funcInd)
            dvct = 0
            for DV in self.DVdict:
                dvname = DV["name"]
                if (DV["type"] == "shape"): 
                    self.shapeGrad[funcInd, dvct] += self.fun3dAim.dynout[funcKey].deriv(dvname)
                    dvct += 1
            
        #update status
        self.cwrite("finished aero mesh contribution to shape DVs\n")


#subclass for tacs steady interface
class TACSinterface(TacsSteadyInterface):
    def __init__(self, comm, tacs_comm, model, n_tacs_procs, datFile, structDVs):
        super(TACSinterface,self).__init__(comm, tacs_comm, model)

        assembler = None
        nodes = None
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

            tInput = structDVs

            def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
                #we need to investigate ordering here, because propID does not match this
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


            #Tim's fix, based on tacs/pyMeshloader.py
            #print("n_tacs_procs={}\n".format(n_tacs_procs),flush=True)
            if (n_tacs_procs > 1):
                ownerRange = assembler.getOwnerRange()
                minNodes = ownerRange[tacs_comm.rank]+1
                maxNodes = ownerRange[tacs_comm.rank+1]+1
                nodes = np.arange(minNodes, maxNodes, dtype=int)
                #print("TACS proc rank {} with nodes={},{}".format(tacs_comm.Get_rank(),minNodes,maxNodes),flush=True)
            else:
                origNodePositions = FEASolver.getOrigNodes()
                nnodes = origNodePositions.shape[0]//3
                nodes = np.arange(1,nnodes+1, dtype=int)

        self._initialize_variables(assembler, struct_id=nodes, thermal_index=6) #
        self.initialize(model.scenarios[0],model.bodies)

    def post_export_f5(self):
        flag = (TACS.OUTPUT_CONNECTIVITY |
                TACS.OUTPUT_NODES |
                TACS.OUTPUT_DISPLACEMENTS |
                TACS.OUTPUT_STRAINS |
                TACS.OUTPUT_STRESSES |
                TACS.OUTPUT_EXTRAS)
        f5 = TACS.ToFH5(self.assembler, TACS.BEAM_OR_SHELL_ELEMENT, flag)
        tacsOutFile = os.path.join(os.getcwd(), "funtofem","TACSoutput.f5")
        f5.writeToFile(tacsOutFile)

##----------Supporting Methods and Classes --------------##

def readnprocs(root_dir=None):
    #get root dir
    if (root_dir is None):
        root_dir = os.getcwd()

    #read the nprocs from run.pbs
    runpbs = os.path.join(root_dir, "run.pbs")
    hdl = open(runpbs, "r")
    lines = hdl.readlines()
    hdl.close()

    #line with nprocs has this format:
    #PBS -l select=4:ncpus=48:mpiprocs=48
    #assume last two numbers are always the same

    ntotprocs = 0
    for line in lines:
        if ("#PBS -l select=" in line):
            chunks = line.split("=")
            # 4:ncpus
            chunk1 = chunks[1]
            subchunk1 = chunk1.split(":")[0]
            ncpus = int(subchunk1)
            # 48:mpiprocs
            chunk2 = chunks[2]
            subchunk2 = chunk2.split(":")[0]
            nprocs = int(subchunk2)

            #total number of procs is product of these two
            ntotprocs = ncpus * nprocs

    return ntotprocs    


def writeInput(DVdict, functions, mode="adjoint", eps=None, x_direction=None):
    # script to write the design variables to funtofem
    # sends them in the file funtofem.in
    # also tells which kind of analysis to perform, etc.

    #make sure funtofem folder exists
    funtofemFolder = os.path.join(os.getcwd(), "funtofem")
    if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

    #make funtofem input file
    inputFile = os.path.join(funtofemFolder, "funtofem.in")
    inputHandle =  open(inputFile, "w")

    #write the real or complex mode, aka adjoint or complex_step
    line = "mode,{}\n".format(mode)
    inputHandle.write(line)

    #write the analysis functions to be requested
    for function in functions:
        line = "function,{}\n".format(function)
        inputHandle.write(line)

    #overwrite values in dict with the input x
    names = []
    values = []

    nDV = 0
    for DV in DVdict:
        #get attributes
        name = DV["name"]
        dvType = DV["type"]
        capsGroup = DV["capsGroup"]

        #get value
        if (DV["active"]):
            nDV += 1

            value = DV["value"]
            #value = round(value, 5)

            #append values and names
            names.append(name)
            values.append(value)

            #make the line strings
            if (dvType == "shape"):
                line = "{},{},{:5f}\n".format(name, dvType, value)
            elif (dvType == "struct"):
                line = "{},{},{},{:5f}\n".format(name, dvType, capsGroup, value)   

            #write the line to the input file
            inputHandle.write(line)


    #if mode is complex_step then also write the direction of complex step
    if (mode == "complex_step"):
        #if epsilon and direction are unspecified make default values
        if (eps is None): eps = 1.0e-30
        if (x_direction is None):
            x_direction = np.random.rand(nDV)
            x_direction = x_direction / np.linalg.norm(x_direction)

        #write the epsilon
        line = "epsilon,{}\n".format(eps)
        inputHandle.write(line)

        #write the direction vector
        inputHandle.write("x_dir,")
        for i in range(len(x_direction)):
            inputHandle.write("{}".format(x_direction[i]))
            if (i < len(x_direction)-1):
                inputHandle.write(",")

    #close the input file
    inputHandle.close()     

def readOutput(DVdict,mode="adjoint"):
    
    #read functions, gradients from funtofem call
    #make sure funtofem folder exists
    funtofemFolder = os.path.join(os.getcwd(), "funtofem")
    if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

    #make funtofem input file
    outputFile = os.path.join(funtofemFolder, "funtofem.out")
    if (os.path.exists(outputFile)):
        #if the output file was written, the analysis ran successfully
        success = True

        #so read in the outputs
        outputHandle =  open(outputFile, "r")

        lines = outputHandle.readlines()
            
        if (mode == "adjoint"):

            functions = {}
            gradients = {}

            #read function values and gradients
            ifunc = -1
            firstLine = True
            for line in lines:
                if ("func" in line):
                    #func, name, value
                    ifunc += 1
                    parts = line.split(",")
                    functionName = parts[1]
                    value = float(parts[2])

                    functions[functionName] = value
                    gradients[functionName] = np.zeros((nDV))
                    iDV = 0

                elif ("grad" in line):
                    #grad,dvname,deriv_i
                    parts = line.split(",")
                    dvname = parts[1]
                    deriv = float(parts[2])

                    #find the DVind of that design variable (assuming out of order)
                    for DV in DVdict:
                        if (DV["name"] == dvname): ind = DV["opt_ind"]
                    #store the gradient at that index
                    gradients[functionName][ind] = deriv

                elif (firstLine):
                    firstLine = False
                    #it's the first line
                    parts = line.split(",")
                    nfunc = int(parts[0])
                    nDV = int(parts[1])

            #close the file
            outputHandle.close()

            #return the functions and gradients
            return functions, gradients

        elif (mode == "forward"):
            functions = {}

            #read function values only
            ifunc = -1
            firstLine = True
            for line in lines:
                if ("func" in line):
                    #func, name, value
                    ifunc += 1
                    parts = line.split(",")
                    functionName = parts[1]
                    value = float(parts[2])

                    functions[functionName] = value
                    iDV = 0

                elif (firstLine):
                    firstLine = False
                    #it's the first line
                    parts = line.split(",")
                    nfunc = int(parts[0])
                    nDV = int(parts[1])

            #close the file
            outputHandle.close()

            #return the functions and gradients
            return functions

        elif (mode == "complex_step"):
            
            #read the funtofem.out with the following format
            #nfunc
            #func,fname,real,value,imag,value
            functions = {}

            for line in lines:
                chunks = line.split(",")
                if (len(chunks) == 1):
                    nfunc = chunks[0].strip()
                elif ("func" in line):
                    name = chunks[1]
                    real_part = float(chunks[3])
                    imag_part = float(chunks[5].strip())

                    functions[name] = real_part + imag_part * 1j

            #return the function values
            return functions


    else:
        sys.exit("Error: Funtofem Analysis failed, no funtofem.out file was created\n")


class Optimize():
    #class to run pyoptsparse optimize on the outside of funtofem

    def __init__(self, DVdict, optimizationMode):
        #set the DV dict here
        self.DVdict = DVdict

        #option: "structural", "full"
        self.optimizationMode = optimizationMode

        #root directory
        self.root_dir = os.getcwd()

        #caps2fun directory
        self.caps2fun_dir = os.environ["CAPS2FUN"]

        #other information such as analysis type, mesher type is  setup in funtofem.py
        self.n_procs = readnprocs()

        #make status file
        optimization_folder = os.path.join(self.root_dir, "optimization")
        if (not(os.path.exists(optimization_folder))): os.mkdir(optimization_folder)
        statusFile = os.path.join(optimization_folder, "opt_status.out")
        self.status = open(statusFile, "w")

        #iterations
        self.iteration = 1

        self.maxStress = None

        #used DVs

        self.cwrite("Aerothermoelastic Optimization with FuntoFem and ESP/CAPS\n")
        self.cwrite("\tDesign Problem: NACA Symmetric Wing\n")
        self.cwrite("Authors: Sean Engelstad, Brian Burke, Sejal Sahu, Graeme Kennedy\n")
        self.cwrite("\tGeorgia Tech SMDO Lab April 2022\n")
        self.cwrite("----------------------------")
        self.cwrite("----------------------------\n")

    def cwrite(self, msg):
        self.status.write(msg)
        self.status.flush()

    def roundVec(self, vec):
        for i in range(len(vec)):
            vec[i] = round(vec[i], 5)
        return vec

    def writeInput(self, x):
        #update design variable dict
        self.updateDVdict(x)

        if (self.optimizationMode == "structural"):
            funcs = ["ksfailure", "mass"]
        elif (self.optimizationMode == "full"):
            funcs = ["ksfailure","cl","cd","mass"]

        writeInput(self.DVdict, funcs)

        #update status to the inputs that were run
        self.cwrite("----------------------------")
        self.cwrite("----------------------------\n")
        self.cwrite("Global Iteration #{}\n".format(self.iteration))
        self.cwrite("\tDV names {}\n".format(self.names))
        self.cwrite("\tDV values {}\n".format(self.values))              

    def updateDVdict(self,x):
        tempDict = []
        self.names = []
        self.values = []
        for DV in self.DVdict:

            #overwrite the value
            if (DV["active"]):
                name = DV["name"]
                DV["value"] = float(x[name])
                
                self.names.append(DV["name"])
                self.values.append(DV["value"])
            tempDict.append(DV)

        print(tempDict)
        
        #overwrite DVdict
        self.DVdict = tempDict

    def callFuntoFem(self):
        #call funtofem via a system call (os.system)
        #the reason for this is fun3d can't be run twice in the system python script
        #this gets around that issue

        #update status
        self.cwrite("\tRunning F2F... ")

        #instead of bash, do system call inside of this python script
        callMessage = "mpiexec_mpt -n {} python $CAPS2FUN/caps2fun/caps2fun.py 2>&1 > ./funtofem/output.txt".format(self.n_procs)
        os.system(callMessage)

        #update status that F2F finished
        self.cwrite("finished F2F!\n")
        #struct DV settings

    def readOutput(self):
        #read functions, gradients from funtofem call
        #make sure funtofem folder exists
        funtofemFolder = os.path.join(os.getcwd(), "funtofem")
        if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

        #make funtofem input file
        self.outputFile = os.path.join(funtofemFolder, "funtofem.out")
        if (os.path.exists(self.outputFile)):
            #if the output file was written, the analysis ran successfully
            self.success = True

            #so read in the outputs
            outputHandle =  open(self.outputFile, "r")

            lines = outputHandle.readlines()
            

            self.functions = {}
            self.gradients = {}

            #read function values and gradients
            ifunc = -1
            firstLine = True
            for line in lines:
                if ("func" in line):
                    #func, name, value
                    ifunc += 1
                    parts = line.split(",")
                    functionName = parts[1]
                    value = float(parts[2])

                    self.functions[functionName] = value
                    self.gradients[functionName] = np.zeros((nDV))
                    iDV = 0

                elif ("grad" in line):
                    #grad,dvname,deriv_i
                    parts = line.split(",")
                    dvname = parts[1]
                    deriv = float(parts[2])

                    #find the DVind of that design variable (assuming out of order)
                    for DV in self.DVdict:
                        if (DV["name"] == dvname):
                            ind = DV["opt_ind"]
                    self.gradients[functionName][ind] = deriv
                    iDV += 1

                elif (firstLine):
                    firstLine = False
                    #it's the first line
                    parts = line.split(",")
                    nfunc = int(parts[0])
                    nDV = int(parts[1])

            #close the file
            outputHandle.close()

            #store the output file in optimization folder
            optimization_folder = os.path.join(self.root_dir, "optimization")
            if (not(os.path.exists(optimization_folder))): os.mkdir(optimization_folder)
            iteration_folder = os.path.join(optimization_folder, "iteration" + str(self.iteration))
            if (not(os.path.exists(iteration_folder))): os.mkdir(iteration_folder)
            for filename in ["funtofem.in", "funtofem.out", "TACSoutput.f5"]:
                src = os.path.join(self.root_dir, "funtofem", filename)
                chunks = filename.split(".")
                newfilename = chunks[0] + str(self.iteration) + "." + chunks[1]
                dest = os.path.join(iteration_folder, newfilename)
                shutil.copy(src, dest)

            #delete output file
            os.remove(self.outputFile)

            #update status
            self.cwrite("\tRead output files for func, sens\n")

            #since didn't fail update iteration count
            self.iteration += 1
            self.fail = 0

        else:
            #otherwise if the output file was not written, it failed
            #possibly due to negative volume or something
            #can try and change parameters here or pass information to change settings and rerun
            self.success = False
            self.cwrite("\tAnalysis Failed\n")
            self.fail = 1
        

        return self.fail

    def deleteF2Ffiles(self):
        os.remove(self.inputFile)
        os.remove(self.outputFile)

    def objCon(self, x):
        #py opt sparse function evaluator

        #write input, call funtofem, and get results
        self.writeInput(x)
        self.callFuntoFem()
        self.readOutput()

        #self.deleteF2Ffiles()

        #objective, constraint functions
        funcs = {}

        if (self.optimizationMode == "structural"):
            funcs["obj"] = self.functions["mass"]
            funcs["con"] = self.functions["ksfailure"]
            objname = "mass"
            conname = "ksfailure"
        elif (self.optimizationMode == "full"):
            #TBD some kind of min fuel burn thing
            funcs["obj"] = 1
            funcs["con"] = 1
            objname = ""
            conname = ""

        #write objective function
        self.cwrite("\t{} Obj = {}\n".format(objname,funcs["obj"]))
        self.cwrite("\t{} Con = {}\n".format(conname,funcs["con"]))

        return funcs, self.fail

    def objGrad(self, x, funcs):
        #pyoptsparse gradient function

        if (self.optimizationMode == "structural"):
            objGrad = self.gradients["mass"]
            conGrad = self.gradients["ksfailure"]
            objname = "mass"
            conname = "ksfailure"
        elif (self.optimizationMode == "full"):
            #TBD some kind of min fuel burn thing
            #objGrad = self.gradients["mass"]
            #conGrad = self.gradients["ksfailure"]
            objname = ""
            conname = ""
        

        #self.cwrite("\t {} Obj grad = {}\n".format(objname,objGrad))
        #self.cwrite("\t{} Con grad = {}\n".format(conname,conGrad))

        sens = {}
        iDV = 0
        objSens = {}
        conSens = {}
        for DV in self.DVdict:
            name = DV["name"]
            if (DV["active"]):
                objSens[name] = objGrad[iDV]
                conSens[name] = conGrad[iDV]
                iDV += 1

        sens["obj"] = objSens
        sens["con"] = conSens

        return sens, self.fail

class Test():
    def __init__(self, DVdict, functions=None):

        self.root_dir = os.getcwd()

        #copy the DV dict
        self.DVdict = DVdict

        #set the num procs
        self.n_procs = readnprocs()

        #set the functions to check
        self.functions = functions
        if (functions is None): 
            self.functions = ["ksfailure","cl","cd","mass"]

        #count the number of active DV
        self.nDV = 0
        for DV in DVdict:
            if (DV["active"]): self.nDV += 1

    def derivativeTest(self):

        #run the adjoint
        self.runAdjoint()

        #run the complex step
        self.runComplexStep()

        #compare the directional derivatives and write to file
        self.writeResults()

    def callCaps2fun(self):
        #call caps2fun through funtofem analysis

        #run funtofem analysis
        callMessage = "mpiexec_mpt -n {} python $CAPS2FUN/caps2fun/caps2fun.py 2>&1 > ./funtofem/output.txt".format(self.n_procs)
        os.system(callMessage)

    def multiForward(self, nruns):
        #run the forward analysis multiple times and store the funtofem.output files in a data folder
        
        funtofemFolder = os.path.join(self.root_dir, "funtofem")
        dataFolder = os.path.join(funtofemFolder, "data")
        if (not(os.path.exists(funtofemFolder))): mkdir(funtofemFolder)
        if (not(os.path.exists(dataFolder))): mkdir(dataFolder)

        for irun in range(nruns):
            
            #call the forward analysis
            self.runForward()

            #copy the funtofem.out file to the data folder
            src = os.path.join(funtofemFolder, "funtofem.out")
            filename = "funtofem" + str(irun+1) + ".out"
            dest = os.path.join(dataFolder, filename)
            shutil.copy(src, dest)
            

    def runForward(self):
        
        mode = "forward"

        #write the F2F input file for adjoint mode
        writeInput(self.DVdict, self.functions, mode=mode)

        #turnoff complex mode, this prob doesn't work
        #os.system("export CMPLX_MODE=0")

        #run funtofem
        self.callCaps2fun()

        self.funcs = readOutput(self.DVdict, mode=mode)

        return self.funcs

    def runAdjoint(self):
        #run the adjoint mode

        #write the F2F input file for adjoint mode
        writeInput(self.DVdict, self.functions, mode="adjoint")

        #turnoff complex mode, this prob doesn't work
        #os.system("export CMPLX_MODE=0")

        #run funtofem
        self.callCaps2fun()

        #read the output file
        self.adjoint_funcs, self.adjoint_grads = readOutput(self.DVdict, mode="adjoint")

    def runComplexStep(self,h=1e-30):
        #run funtofem in complex step mode

        #generate random perturbation for complex_step check
        x_dir = np.random.rand(self.nDV)
        x_dir = x_dir / np.linalg.norm(x_dir)

        #write the F2F input file for complex mode
        writeInput(self.DVdict, self.functions, mode="complex_step", eps=h, x_direction=x_dir)

        #run funtofem
        self.callCaps2fun()

        #read the output file in complex mode, complex_step.out
        self.complex_funcs = readOutput(self.DVdict,mode="complex_step")

    def writeResults(self):
        #compare directional derivatives of adjoint vs complex step

        #make a derivative_check.out file in funtofem folder
        funtofemFolder = os.path.join(os.getcwd(), "funtofem")
        deriv_file = os.path.join(funtofemFolder, "derivative_check.out")

        deriv_handle = open(deriv_file, "w")

        fileExists = os.path.exists(deriv_file)

        #write a an introductory line to the file
        if (not(fileExists)):
            deriv_handle.write("Funtofem Derivative Check, Adjoint vs Complex Step\n")
            deriv_handle.write("----------------------------")
            deriv_handle.write("----------------------------\n")

        adjoint_dderiv = {}
        complex_dderiv = {}

        #for each function in the analysis
        for function in self.functions:
            
            #write the function name "func,name"
            deriv_handle.write("func,{}\n".format(function))

            #write the adjoint and complex_step functions (real part) to the file
            deriv_handle.write("\tadjoint func      = {}\n".format(self.adjoint_funcs[function]))
            deriv_handle.write("\tcomplex_step func = {}\n".format(self.complex_funcs[function].real))

            #initialize directional derivative, and get this adjoint_grad
            adjoint_dderiv[function] = 0
            adjoint_grad = self.adjoint_grads[function]

            #loop over each design variable and dot product the perturbation and gradient to get directional derivative
            for i in range(self.nDV):
                adjoint_dderiv[function] += x_dir[i] * adjoint_grad[i]

            #write the adjoint directional derivative to the file
            line = "\tadjoint dderiv     = {}\n".format(adjoint_dderiv[function])
            deriv_handle.write(line)

            #compute complex step derivative, which is also for that direction
            #f(x+hp*1j)/h is the complex step
            complex_dderiv[function] = self.complex_funcs[function].imag

            #write the complex step directional derivative to the file
            line = "\tcomplex_step dderiv = {}\n".format(complex_dderiv[function])
            deriv_handle.write(line)

            #compute the absolute error
            absoluteError = abs( complex_dderiv[function] - adjoint_dderiv[function] )

            #write the absolute error to the file
            line = "\tabsolute error      = {}\n".format(absoluteError)
            deriv_handle.write(line)

            #compute the relative error
            if (abs(complex_dderiv[function]) < 1.0e-15):
                relativeError = None
            else:
                relativeError = absoluteError / abs(complex_dderiv[function])

            #write relative error to a file
            line = "\relative error       = {}\n".format(relativeError)
            deriv_handle.write(line)


        #close the derivative_check.out file
        deriv_handle.close()

##----------Outside of class, run cases------------------##

if (inputExists):
    #run Caps2fun
    myrun = Caps2Fun()
    myrun.runF2F()