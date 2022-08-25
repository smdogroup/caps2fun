import caps2tacs
import numpy as np
import os

panel_problem = caps2tacs.CapsStruct.default(csmFile="nacaWing.csm")
tacs_aim = panel_problem.tacsAim
egads_aim = panel_problem.egadsAim

aluminum = caps2tacs.Isotropic.aluminum()
tacs_aim.add_material(material=aluminum)

# add wing root constraint
constraint = caps2tacs.ZeroConstraint(name="fixRoot", caps_constraint="wingRoot")
tacs_aim.add_constraint(constraint=constraint)

load = caps2tacs.GridForce(name="load1", caps_load="fullWing", direction=[0.,0.,1.], magnitude=1.0E3)    
tacs_aim.add_load(load=load)

thick_idx = 0
nribs = 16
for rib_idx in range(1,nribs+1):
    thick_DV = caps2tacs.ThicknessVariable(name=f"thick{thick_idx}", caps_group=f"rib{rib_idx}", value=0.2-0.005*rib_idx, material=aluminum)
    tacs_aim.add_variable(variable=thick_DV)
    thick_idx += 1

nspars = 2
for spar_idx in range(1,nspars+1):
    thick_DV = caps2tacs.ThicknessVariable(name=f"thick{thick_idx}", caps_group=f"spar{spar_idx}", value=0.4-0.08*spar_idx, material=aluminum)
    tacs_aim.add_variable(variable=thick_DV)
    thick_idx += 1

nOML = nribs-1
for OML_idx in range(1,nOML+1):
    thick_DV = caps2tacs.ThicknessVariable(name=f"thick{thick_idx}", caps_group=f"OML{OML_idx}", value=0.05-0.002*OML_idx, material=aluminum)
    tacs_aim.add_variable(variable=thick_DV)
    thick_idx += 1

tacs_aim.setup_aim()
egads_aim.set_mesh(edge_pt_min=5, edge_pt_max=10, global_mesh_size=0.50, max_surf_offset=0.01, max_dihedral_angle=15)

# make a pytacs function
pytacs_function = caps2tacs.MassStressTransient(t0=0.0, tf=10.0, num_steps=10, amplitude=lambda t : np.sin(0.5*t))

caps_tacs = caps2tacs.CapsTacs(
    name="naca_wing_struct", tacs_aim=tacs_aim, 
    egads_aim=egads_aim, pytacs_function=pytacs_function, write_solution=True, view_plots=False
    )
caps_tacs.analysis()

# convert all f5 to vtk in analysis dir 
for filename in os.listdir(tacs_aim.analysis_dir):
    filepath = os.path.join(tacs_aim.analysis_dir, filename)
    os.system(f"~/git/tacs/extern/f5tovtk/f5tovtk {filepath}")