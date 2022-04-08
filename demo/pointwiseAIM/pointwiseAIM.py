#pointwise aim demo

import pyCAPS
import os

def run_pointwise(pointwise):
    #run AIM pre-analysis
    pointwise.preAnalysis()

    #move to test directory
    currentDir = os.getcwd()
    os.chdir(pointwise.analysisDir)

    CAPS_GLYPH = os.environ["CAPS_GLYPH"]
    for i in range(60):
        os.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")

        if os.path.isfile('caps.GeomToMesh.gma') and os.path.isfile('caps.GeomToMesh.ugrid'): break
    
    os.chdir(currentDir)

    #run AIM postanalysis
    pointwise.postAnalysis()

#------------------------------------------------#

#load the csm file
filename = os.path.join("naca_small.csm")
myProblem = pyCAPS.Problem(problemName = "workDir", capsFile = filename, outLevel = 1)

#create pointwise aim
pointwise = myProblem.analysis.create(aim = "pointwiseAIM", name = "pointwise")

#Dump VTK files for visualization
pointwise.input.Proj_Name = "Transport"
pointwise.input.Mesh_Format = "VTK"

#Connector Level
pointwise.input.Connector_Turn_Angle = 10
pointwise.input.Connector_Turn_Angle_Hard = 70
pointwise.input.Connector_Source_Spacing = True

#Domain level
pointwise.input.Domain_Algorithm = "AdvancingFront"
pointwise.input.Domain_Max_Layers = 15
pointwise.input.Domain_TRex_ARLimit = 40.0
pointwise.input.Domain_Decay = 0.8

#Block level
pointwise.input.Block_Boundary_Decay = 0.8
pointwise.input.Block_Edge_Max_Growth_Rate = 1.5

#Execute pointwise
run_pointwise(pointwise)

#view the surface tesselation
pointwise.geometry.view()
