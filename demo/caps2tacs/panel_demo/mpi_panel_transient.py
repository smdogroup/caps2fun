import caps2tacs
import os
import numpy as np
from mpi4py import MPI

# all the capsWrapper package modules/classes are also available 
# in caps2tacs since caps2tacs package imports the whole capsWrapper package\

comm = MPI.COMM_WORLD

tacs_aim = None
egads_aim = None

if comm.Get_rank() == 0:

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

    tacs_aim.setup_aim()
    egads_aim.set_mesh()

pytacs_function = caps2tacs.MassStressTransient(t0=0.0, tf=10.0, num_steps=100, amplitude=lambda t : np.sin(0.5*t))

# start a caps tacs main problem
caps_tacs = caps2tacs.CapsTacs(name="transient_panel", 
tacs_aim=tacs_aim, egads_aim=egads_aim, pytacs_function=pytacs_function, comm=comm,
compute_gradients=False, view_plots=True)

caps_tacs.analysis()


