import capsManager
from tacs.pytacs import pyTACS
from tacs import functions
import os

panel_problem = capsManager.CapsStruct.default(csmFile="nacaWing.csm")
tacs_aim = panel_problem.tacsAim
egads_aim = panel_problem.egadsAim

OML_madeupium = capsManager.Isotropic.madeupium()
OML_madeupium.name = "OML-madeupium"
tacs_aim.add_material(material=OML_madeupium)

rib_madeupium = capsManager.Isotropic.madeupium()
rib_madeupium.name = "rib-madeupium"
rib_madeupium.young_modulus *= 1.5
tacs_aim.add_material(material=rib_madeupium)


spar_madeupium = capsManager.Isotropic.madeupium()
spar_madeupium.name = "spar-madeupium"
spar_madeupium.young_modulus *= 2.0
tacs_aim.add_material(material=spar_madeupium)

# add wing root constraint
constraint = capsManager.ZeroConstraint(name="fixRoot", caps_constraint="wingRoot")
tacs_aim.add_constraint(constraint=constraint)

load = capsManager.GridForce(name="load1", caps_load="fullWing", direction=[1.,0.,0.], magnitude=1.0E3)    
tacs_aim.add_load(load=load)

thick_idx = 0
nribs = 16
for rib_idx in range(1,nribs+1):
    thick_DV = capsManager.ThicknessVariable(name=f"thick{thick_idx}", caps_group=f"rib{rib_idx}", value=0.2-0.005*rib_idx, material=rib_madeupium)
    tacs_aim.add_variable(variable=thick_DV)
    thick_idx += 1

nspars = 2
for spar_idx in range(1,nspars+1):
    thick_DV = capsManager.ThicknessVariable(name=f"thick{thick_idx}", caps_group=f"spar{spar_idx}", value=0.4-0.08*spar_idx, material=spar_madeupium)
    tacs_aim.add_variable(variable=thick_DV)
    thick_idx += 1

nOML = nribs-1
for OML_idx in range(1,nOML+1):
    thick_DV = capsManager.ThicknessVariable(name=f"thick{thick_idx}", caps_group=f"OML{OML_idx}", value=0.05-0.002*OML_idx, material=OML_madeupium)
    tacs_aim.add_variable(variable=thick_DV)
    thick_idx += 1

despmtrs = panel_problem.geometry.despmtr.keys()
print(despmtrs)
for despmtr in despmtrs:
    shape_var = capsManager.ShapeVariable(name=despmtr)
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





