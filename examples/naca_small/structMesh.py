##-----------------import stuff--------------------------##
import os
import pyCAPS
from mpi4py import MPI
from tacs.pytacs import pyTACS
from tacs import functions

##-------------------------------------------------------##
filename = os.path.join("wing.egads")
caps = pyCAPS.Problem(problemName = "myCAPS",
                    capsFile = filename,
                    outLevel = 1)

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
# Solve each structural problem and write solutions
func = {}; sens = {}
for caseID in SPs:
    SPs[caseID].solve()
    print("finished pytacs solve")
    SPs[caseID].evalFunctions(func,evalFuncs=evalFuncs)
    #print("finished pytacs funcs")
    SPs[caseID].evalFunctionsSens(sens,evalFuncs=evalFuncs)
    #print("finished pytacs sens")
    SPs[caseID].writeSolution(outputDir=tacsAim.analysisDir)
    #print("finished pytacs file")

curDir = os.getcwd()
os.chdir(tacsAim.analysisDir)
os.system("f5tovtk *.f5")
os.chdir(curDir)