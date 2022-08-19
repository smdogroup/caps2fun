import caps2tacs
from tacs.pytacs import pyTACS
from tacs import functions
import os, sys

# all the capsWrapper package modules/classes are also available 
# in caps2tacs since caps2tacs package imports the whole capsWrapper package

problem = caps2tacs.CapsStruct.default(csmFile="purdue_p.csm")
tacs_aim = problem.tacsAim
egads_aim = problem.egadsAim

madeupium = caps2tacs.Isotropic.madeupium()
tacs_aim.add_material(material=madeupium)

constraint = caps2tacs.ZeroConstraint(name="bottom", caps_constraint="bottom")
tacs_aim.add_constraint(constraint=constraint)

for group in ["front", "back", "boundary"]:
    thick_DV = caps2tacs.ThicknessVariable(name=group, caps_group=group, value=0.01, material=madeupium)
    tacs_aim.add_variable(variable=thick_DV)

load = caps2tacs.GridForce(name="load1", caps_load="front", direction=[0,0,-1.0], magnitude=1.0E3)    
tacs_aim.add_load(load=load)

for despmtr in problem.design_parameters:
    shape_var = caps2tacs.ShapeVariable(name=despmtr)
    tacs_aim.add_variable(variable=shape_var)

tacs_aim.setup_aim()
egads_aim.set_mesh(edge_pt_min=40, edge_pt_max=50, tess_params=[0.01,0.01,15])

# start a caps tacs main problem
caps_tacs = caps2tacs.CapsTacs(tacs_aim=tacs_aim, egads_aim=egads_aim)

#initialize pytacs with that data file
FEASolver = pyTACS(caps_tacs.dat_file)
    
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





