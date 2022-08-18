import capsWrapper
from tacs.pytacs import pyTACS
from tacs import functions
import os

panel_problem = capsWrapper.CapsStruct.default(csmFile="panel.csm")
tacs_aim = panel_problem.tacsAim
egads_aim = panel_problem.egadsAim

madeupium = capsWrapper.Isotropic.madeupium()
tacs_aim.add_material(material=madeupium)


constraint = capsWrapper.ZeroConstraint(name="fixEdge", caps_constraint="edge")
tacs_aim.add_constraint(constraint=constraint)

load = capsWrapper.GridForce(name="load1", caps_load="plate", direction=[0,0,-1.], magnitude=1.0E3)    
tacs_aim.add_load(load=load)

thick_DV = capsWrapper.ThicknessVariable(name="thick", caps_group="plate", value=0.01, material=madeupium)
tacs_aim.add_variable(variable=thick_DV)

despmtrs = ["plateLength", "plateWidth"]
for despmtr in despmtrs:
    shape_var = capsWrapper.ShapeVariable(name=despmtr)
    tacs_aim.add_variable(variable=shape_var)

tacs_aim.setup_aim()

built_tacs_aim = tacs_aim.aim
built_egads_aim = egads_aim.aim

built_tacs_aim.input["Mesh"].link(built_egads_aim.output["Surface_Mesh"])

# run the tacs aim preanalysis
built_tacs_aim.preAnalysis()

dat_file = os.path.join(built_tacs_aim.analysisDir, built_tacs_aim.input.Proj_Name + '.dat')

#initialize pytacs with that data file
FEASolver = pyTACS(dat_file)
    
# Set up TACS Assembler
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

# Solve each structural problem and write solutions
func = {}; sens = {}
for caseID in SPs:
    SPs[caseID].solve()
    SPs[caseID].evalFunctions(func,evalFuncs=evalFuncs)
    SPs[caseID].evalFunctionsSens(sens,evalFuncs=evalFuncs)
    SPs[caseID].writeSolution(outputDir=os.path.dirname(__file__))





