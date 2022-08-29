import caps2tacs
import sys
import numpy as np
from mpi4py import MPI

# all the capsManager package modules/classes are also available 
# in caps2tacs since caps2tacs package imports the whole capsManager package

# pyoptsparse documentation: https://mdolab-pyoptsparse.readthedocs-hosted.com/en/latest/

comm = MPI.COMM_WORLD

tacs_aim = None
egads_aim = None

if comm.Get_rank() == 0:

    problem = caps2tacs.CapsStruct.default(csmFile="generic_missile.csm")
    tacs_aim = problem.tacsAim
    egads_aim = problem.egadsAim

    aluminum = caps2tacs.Isotropic.aluminum()
    tacs_aim.add_material(material=aluminum)

    constraint = caps2tacs.ZeroConstraint(name="bottom", caps_constraint="root")
    tacs_aim.add_constraint(constraint=constraint)

    for group in ["MissileOML"]:
        thick_DV = caps2tacs.ThicknessVariable(name=group, caps_group=group, value=0.1, material=aluminum)
        tacs_aim.add_variable(variable=thick_DV)

    load = caps2tacs.GridForce(name="load1", caps_load="fuselage", direction=[0,0,-1.0], magnitude=1.0E2)    
    tacs_aim.add_load(load=load)

    tacs_aim.setup_aim()
    egads_aim.set_mesh(edge_pt_min=30, edge_pt_max=40, global_mesh_size=0.01)

#print([dv.name for dv in tacs_aim.shape_design_variables])

# make a pytacs function
pytacs_function = caps2tacs.MassStressTransient(t0=0.0, tf=10.0, num_steps=100, amplitude=lambda t : np.sin(0.5*t))

caps_tacs = caps2tacs.CapsTacs(
    name="missile_transient", tacs_aim=tacs_aim, 
    egads_aim=egads_aim, comm=comm,
    pytacs_function=pytacs_function, write_solution=True, 
    view_plots=True
    )
caps_tacs.analysis()