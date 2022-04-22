#------------------------------------------------------------------------------#
import pyCAPS
import os

#------------------------------------------------------------------------------#


def run_pointwise(pointwise):
    #run AIM pre-analysis
    pointwise.preAnalysis()

    #move to test directory
    currentDir = os.getcwd()
    os.chdir(pointwise.analysisDir)

    CAPS_GLYPH = os.environ["CAPS_GLYPH"]
    for i in range(1):
        os.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")

        #if os.path.isfile('caps.GeomToMesh.gma') and os.path.isfile('caps.GeomToMesh.ugrid'): break
    
    os.chdir(currentDir)

    #run AIM postanalysis
    pointwise.postAnalysis()

#------------------------------------------------------------------------------#
# Load CSM file
#os.system("serveCSM naca_small.csm")

filename = os.path.join("cfd.egads")
caps = pyCAPS.Problem(problemName = "myCAPS",
                    capsFile = filename,
                    outLevel = 1)

# Create pointwise aim
pointwise = caps.analysis.create(aim = "pointwiseAIM",
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
        "Farfield": {"bcType":"farfield"}}

# Execute pointwise
run_pointwise(pointwise)

# View the surface tessellation
pointwise.geometry.view()
#os.system("pointwise ")

