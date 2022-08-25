import capsManager
from tacs.pytacs import pyTACS
from tacs import functions
import os

panel_problem = capsManager.CapsStruct.default(csmFile="stiffPanel4.csm")
tacs_aim = panel_problem.tacsAim
egads_aim = panel_problem.egadsAim

madeupium = capsManager.Isotropic.madeupium()
tacs_aim.add_material(material=madeupium)

stronger_madeupium = capsManager.Isotropic.madeupium()
stronger_madeupium.name = "stronger-madeupium"
stronger_madeupium.young_modulus *= 5.0
tacs_aim.add_material(material=stronger_madeupium)

# two constraints for this problem
for idx in range(2):
    constraint = capsManager.ZeroConstraint(name=f"fixEdge{idx}", caps_constraint=f"edge{idx+1}")
    tacs_aim.add_constraint(constraint=constraint)

load = capsManager.GridForce(name="load1", caps_load="plate", direction=[0,-1,-1.], magnitude=1.0E3)    
tacs_aim.add_load(load=load)

n_plates = 8
for plate_idx in range(1,n_plates+1):
    thick_DV = capsManager.ThicknessVariable(name=f"thick{plate_idx}", caps_group=f"plate{plate_idx}", value=0.01+abs(0.1-0.03*plate_idx), material=madeupium)
    tacs_aim.add_variable(variable=thick_DV)

n_stiffeners = 4
for stiff_idx in range(1,n_stiffeners+1):
    thick_DV = capsManager.ThicknessVariable(name=f"thick{n_plates+1+stiff_idx}", caps_group=f"stiffener{stiff_idx}", value=0.2, material=stronger_madeupium)
    tacs_aim.add_variable(variable=thick_DV)

despmtrs = ["plateLength", "plateWidth", "stiffHeight"]
for despmtr in despmtrs:
    shape_var = capsManager.ShapeVariable(name=despmtr)
    tacs_aim.add_variable(variable=shape_var)

# setup both aims
tacs_aim.setup_aim()
egads_aim.set_mesh()

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





