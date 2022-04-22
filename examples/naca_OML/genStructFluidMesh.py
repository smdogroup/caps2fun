import os, shutil
import pyCAPS
from mpi4py import MPI
from tacs.pytacs import pyTACS
from tacs import functions

##--------set design parameters------##
#taper = 0.5
isHPC = True

##--------struct mesh----------------##

#start CAPS problem, clears stuct folder
#need naca_OML.csm default to cfdOn=0 for struct
caps = pyCAPS.Problem(problemName = "struct",
                    capsFile = "naca_OML_struct.csm",
                    outLevel = 1)
#caps.geometry.despmtr["taper"].value = taper

#initialize egads Aim
egadsAim = caps.analysis.create(aim="egadsTessAIM")

#initialize tacs Aim
tacsAim = caps.analysis.create(aim = "tacsAIM", name = "tacs")


##-----------PyCAPS setup function-------##
#setup function for naca_small.csm
    
#Egads Aim section, for mesh
egadsAim.input.Edge_Point_Min = 5
egadsAim.input.Edge_Point_Max = 10

egadsAim.input.Mesh_Elements = "Quad"

egadsAim.input.Tess_Params = [.25,.01,15]

#increase the precision in the BDF file
tacsAim.input.File_Format = "Large"
tacsAim.input.Mesh_File_Format = "Large"

# Link the mesh
tacsAim.input["Mesh"].link(egadsAim.output["Surface_Mesh"])

# Set analysis type
tacsAim.input.Analysis_Type = "Static"

#materials section    
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3,
                "tensionAllow" :  20.0e7}

tacsAim.input.Material = {"madeupium": madeupium}

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

tacsAim.input.Property = {"rib": ribshell,
"spar" : sparshell,
"OML" : OMLshell}

# constraint section
constraint1 = {"groupName" : "wingRoot",
                "dofConstraint" : 123456}

tacsAim.input.Constraint = {"fixRoot": constraint1}

# Set load
liftload = {"groupName"         : "bottomWing",
    "loadType"          : "GridForce",
    "forceScaleFactor"  : 1.0e2,
    "directionVector"   : [0.0, 1.0, 0.0]}

# Set loads
tacsAim.input.Load = {"lift": liftload }

#return the capsGroups you want to have thickDVs in that order
#[thick1, thick2, thick3]
capsDVgroups = ["rib", "spar", "OML"]

##-----------setup DVs and DVRs-----------##
#setup the Design_Variable and Design_Variable_Relations

#list of design variables, with thick1-thickN for thickness DVs
desvars = ["area","aspect","taper","ctwist","lesweep","dihedral","thick1", "thick2", "thick3"]
nvar = len(desvars)

#where the thickDVs have thick1 is for "rib", thick2 for "spar" etc.
capsGroups = ["rib","spar","OML"]
thickness = [0.01, 0.02, 0.03]

def makeThicknessDV(capsGroup, thickness):
    #thick DV dictionary for Design_Variable Dict
    desvar    = {"groupName" : capsGroup,
            "initialValue" : thickness,
            "lowerBound" : thickness*0.5,
            "upperBound" : thickness*1.5,
            "maxDelta"   : thickness*0.1}
    return desvar
    
def makeThicknessDVR(DVname):
    #thick DV dictionary for Design_Variable_Relation Dict
    DVR = {"variableType": "Property",
    "fieldName" : "T",
    "constantCoeff" : 0.0,
    "groupName" : DVname,
    "linearCoeff" : 1.0}
    return DVR

#make initial DV and DVR dict
DVdict = {}
DVRdict = {}

#add thickDVs and geomDVs to caps
thickCt = 0
for desvar in desvars:
    if ("thick" in desvar):
        #add thickDV entry into DV_Relations and DV Dicts
        DVRdict[desvar] = makeThicknessDVR(desvar)
        DVdict[desvar] = makeThicknessDV(capsGroups[thickCt],thickness[thickCt])
        thickCt += 1
    else: #geomDV, add empty entry into DV dicts
        DVdict[desvar] = {}

#input DVdict and DVRdict into tacsAim
tacsAim.input.Design_Variable = DVdict
tacsAim.input.Design_Variable_Relation = DVRdict

##---------PYTACS to run TACS-------------##
#first run the preanalysis to prepare to run tacs through pytacs
tacsAim.preAnalysis()

#setup MPI COMM for pytacs
comm = MPI.COMM_WORLD

#data file
datFile = os.path.join(tacsAim.analysisDir, tacsAim.input.Proj_Name + '.dat')

#initialize pytacs with that data file
FEASolver = pyTACS(datFile)
    
# Set up TACS Assemblerrom tacs
FEASolver.initialize()

#choose the functions to evaluate
evalFuncs = ['wing_mass', 'ks_vmfailure']

#read the bdf & dat file into pytacs FEAsolver
#SPs represents "StructuralProblems"
SPs = FEASolver.createTACSProbsFromBDF()

# Read in forces from BDF and create tacs struct problems
for caseID in SPs:
    SPs[caseID].addFunction('wing_mass', functions.StructuralMass)
    SPs[caseID].addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5, ksWeight=1000.0)
    #SPs[caseID].addFunction('compliance', functions.Compliance)

curDir = os.getcwd()
meshDir = os.path.join(curDir, "mesh")

# Solve each structural problem and write solutions
func = {}; sens = {}
for caseID in SPs:
    SPs[caseID].solve()
    print("finished pytacs solve")
    SPs[caseID].evalFunctions(func,evalFuncs=evalFuncs)
    SPs[caseID].evalFunctionsSens(sens,evalFuncs=evalFuncs)
    SPs[caseID].writeSolution(outputDir=tacsAim.analysisDir)
    SPs[caseID].writeSolution(outputDir=meshDir)

builtStruct = True


##-----------------Fluid Mesh-------------------------##
#load pointwise if on HPC
#module load pointwise/18.5R1

#pointwise or tetgen style
mesh_style = "tetgen"

if (mesh_style == "pointwise"):
    #save cfd egads version for fluid mesh to read next

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

    #-------------------Start fluid capsProblem------------------------------------------#

    caps2 = pyCAPS.Problem(problemName = "fluid",
                        capsFile = "naca_OML_fluid.csm",
                        outLevel = 1)

    #update despmtr
    #caps2.geometry.despmtr["taper"].value = taper

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
            "Farfield": {"bcType":"Farfield"},
            "Symmetry" : {"bcType" : "SymmetryY"}}

    # Execute pointwise
    run_pointwise(pointwise)

    #ran fluid mesh boolean
    builtFluid = True

    #copy files from pointwise analysis dir to meshDir
    os.system("cp " + pointwise.analysisDir + "/*.ugrid ./mesh/") 

else:
    capsFluid = pyCAPS.Problem(problemName = "fluid",
                        capsFile = "naca_OML_fluid.csm",
                        outLevel = 1)
    egadsFluidAim = capsFluid.analysis.create(aim = "egadsTessAIM", name = "egadsTess")
    tetgenAim = capsFluid.analysis.create(aim = "tetgenAIM", name = "tetgen")
    fun3dAim = capsFluid.analysis.create(aim = "fun3dAIM",
                                    name = "fun3d")

    egadsFluidAim.input.Tess_Params = [1.5, 0.001, 0.5]
    egadsFluidAim.input.Edge_Point_Min = 5
    egadsFluidAim.input.Edge_Point_Max = 30

    #tetgen AIM mesh params
    tetgenAim.input.Preserve_Surf_Mesh = True
    tetgenAim.input["Surface_Mesh"].link(egadsFluidAim.output["Surface_Mesh"])
    tetgenAim.input.Mesh_Format = "AFLR3"

    '''
    meshSizeDict = {}
    meshSizeDict["OML"] = {"tessParams" : [0.25, 0.01, 10.0]}
    meshSizeDict["Farfield"] = {"tessParams" : [7.5, 0.01, 10.0]}
    meshSizeDict["Symmetry"] = {"tessParams" : [1, 0.01, 10.0]}
    tetgenAim.input.Mesh_Sizing = meshSizeDict
    '''

    fun3dAim.input["Mesh"].link(tetgenAim.output["Volume_Mesh"])

    fun3dAim.input.Boundary_Condition = {"OML": {"bcType" : "Inviscid"},
                "Farfield": {"bcType":"Farfield"},
                "Symmetry" : "SymmetryZ"}
    fun3dAim.preAnalysis()