import caps2tacs
from mpi4py import MPI

# all the capsWrapper package modules/classes are also available 
# in caps2tacs since caps2tacs package imports the whole capsWrapper package\

panel_problem = caps2tacs.CapsStruct.default(csmFile="panel.csm")
tacs_aim = panel_problem.tacsAim
egads_aim = panel_problem.egadsAim

madeupium = caps2tacs.Isotropic.madeupium()
tacs_aim.add_material(material=madeupium)


constraint = caps2tacs.ZeroConstraint(name="fixEdge", caps_constraint="edge")
tacs_aim.add_constraint(constraint=constraint)

load = caps2tacs.GridForce(name="load1", caps_load="plate", direction=[0,0,-1.], magnitude=1.0E3)    
tacs_aim.add_load(load=load)

thick_DV = caps2tacs.ThicknessVariable(name="thick", caps_group="plate", value=0.01, material=madeupium)
tacs_aim.add_variable(variable=thick_DV)

despmtrs = ["plateLength", "plateWidth"]
for despmtr in despmtrs:
    shape_var = caps2tacs.ShapeVariable(name=despmtr)
    tacs_aim.add_variable(variable=shape_var)

tacs_aim.setup_aim()
egads_aim.set_mesh()

# start a caps tacs main problem
caps_tacs = caps2tacs.CapsTacs(tacs_aim=tacs_aim, egads_aim=egads_aim)

comm = MPI.COMM_WORLD
solver = caps2tacs.TacsSolver(comm=comm, dat_file=caps_tacs.dat_file)
