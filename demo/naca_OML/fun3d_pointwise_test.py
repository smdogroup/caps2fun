import os
import pyCAPS
from mpi4py import MPI
from tacs.pytacs import pyTACS
from tacs import functions

##--------set design parameters------##
taper = 1.0
isHPC = True
##-----------------Fluid Mesh-------------------------##
#save cfd egads version for fluid mesh to read next

def run_pointwise(pointwise):
    #run AIM pre-analysis
    pointwise.preAnalysis()

    #move to test directory
    currentDir = os.getcwd()
    os.chdir(pointwise.analysisDir)

    CAPS_GLYPH = os.environ["CAPS_GLYPH"]
    for i in range(1):
        if (isHPC):
            os.system("vglrun pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")        
        else:
            os.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")

        #if os.path.isfile('caps.GeomToMesh.gma') and os.path.isfile('caps.GeomToMesh.ugrid'): break
    
    os.chdir(currentDir)

    #run AIM postanalysis
    pointwise.postAnalysis()

#-------------------Start fluid capsProblem------------------------------------------#

caps2 = pyCAPS.Problem(problemName = "fluid",
                    capsFile = "naca_OML_fluid.csm",
                    outLevel = 1)

#update despmtr
caps2.geometry.despmtr["taper"].value = taper

#no need to change here since it reads from egads
#caps2.geometry.despmtr["taper"] = taper

# Create pointwise aim
pointwise = caps2.analysis.create(aim = "pointwiseAIM",
                                    name = "pointwise")


#pointwise.geometry.cfgpmtr["cfdOn"].value = 1

# Dump VTK files for visualization
pointwise.input.Proj_Name   = "TransportWing"
pointwise.input.Mesh_Format = "VTK"

# Connector level
pointwise.input.Connector_Turn_Angle       = 10
pointwise.input.Connector_Prox_Growth_Rate = 1.2
pointwise.input.Connector_Source_Spacing   = True

# Domain level
pointwise.input.Domain_Algorithm    = "AdvancingFront"
pointwise.input.Domain_Max_Layers   = 15
pointwise.input.Domain_Growth_Rate  = 1.25
pointwise.input.Domain_TRex_ARLimit = 40.0
pointwise.input.Domain_Decay        = 0.8
pointwise.input.Domain_Iso_Type = "Triangle"

# Block level
pointwise.input.Block_Boundary_Decay       = 0.8
pointwise.input.Block_Collision_Buffer     = 1.0
pointwise.input.Block_Max_Skew_Angle       = 160.0
pointwise.input.Block_Edge_Max_Growth_Rate = 1.5
pointwise.input.Block_Full_Layers          = 1
pointwise.input.Block_Max_Layers           = 100
pointwise.input.Block_TRexType = "TetPyramid"
#T-Rex cell type (TetPyramid, TetPyramidPrismHex, AllAndConvertWallDoms).

# Set wall spacing for capsMesh == leftWing and capsMesh == riteWing
viscousWall  = {"boundaryLayerSpacing" : 0.001}
pointwise.input.Mesh_Sizing = {"OML": viscousWall,
        "Symmetry" : {"bcType" : "SymmetryY"},
        "Farfield": {"bcType":"Farfield"}}

# Execute pointwise
run_pointwise(pointwise)

#ran fluid mesh boolean
builtFluid = True

##---------------Feedback------------##
if (builtFluid): print("Built fluid mesh!")
