#------------------------------------------------------------------------------#

# Import pyCAPS module
import pyCAPS

# Import os and platform
import os
import platform
import time

#------------------------------------------------------------------------------#

def run_pointwise(pointwise):
    # Run AIM pre-analysis
    pointwise.preAnalysis()

    ####### Run pointwise #################
    currentDirectory = os.getcwd() # Get current working directory
    os.chdir(pointwise.analysisDir)    # Move into test directory

    CAPS_GLYPH = os.environ["CAPS_GLYPH"]
    for i in range(1):
        if "Windows" in platform.system():
            PW_HOME = os.environ["PW_HOME"]
            os.system('"' + PW_HOME + '\\win64\\bin\\tclsh.exe ' + CAPS_GLYPH + '\\GeomToMesh.glf" caps.egads capsUserDefaults.glf')
        else:
            os.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")

        time.sleep(1) # let the harddrive breathe
        if os.path.isfile('caps.GeomToMesh.gma') and os.path.isfile('caps.GeomToMesh.ugrid'): break
        time.sleep(20) # wait and try again

    os.chdir(currentDirectory)     # Move back to top directory
    #######################################

    # Run AIM post-analysis
    pointwise.postAnalysis()

#------------------------------------------------------------------------------#

# Load CSM file
filename = os.path.join("wing.egads")
myProblem = pyCAPS.Problem(problemName = "workDir_08_ViscousWing",
                           capsFile = filename,
                           outLevel = 1)

# Create pointwise aim
pointwise = myProblem.analysis.create(aim = "pointwiseAIM",
                                      name = "pointwise")

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

# Block level
pointwise.input.Block_Boundary_Decay       = 0.8
pointwise.input.Block_Collision_Buffer     = 1.0
pointwise.input.Block_Max_Skew_Angle       = 160.0
pointwise.input.Block_Edge_Max_Growth_Rate = 1.5
pointwise.input.Block_Full_Layers          = 1
pointwise.input.Block_Max_Layers           = 100

# Set wall spacing for capsMesh == leftWing and capsMesh == riteWing
viscousWall  = {"boundaryLayerSpacing" : 0.001}
pointwise.input.Mesh_Sizing = {"wing": viscousWall}

# Execute pointwise
run_pointwise(pointwise)

# View the surface tessellation
pointwise.geometry.view()
