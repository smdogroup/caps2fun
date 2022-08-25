import caps2tacs
import os, sys
import numpy as np

# all the capsWrapper package modules/classes are also available 
# in caps2tacs since caps2tacs package imports the whole capsWrapper package

panel_problem = caps2tacs.CapsStruct.default(csmFile="multi_panel.csm")
tacs_aim = panel_problem.tacsAim
egads_aim = panel_problem.egadsAim

madeupium = caps2tacs.Isotropic.madeupium()
tacs_aim.add_material(material=madeupium)

for iplate in range(3):
    group_name = f"plate{iplate+1}"

    if not iplate == 1:
        constraint = caps2tacs.ZeroConstraint(name=f"fixEdge{iplate}", caps_constraint=group_name)
        tacs_aim.add_constraint(constraint=constraint)

    thick_DV = caps2tacs.ThicknessVariable(name="thick", caps_group=group_name, value=0.01+0.005*iplate, material=madeupium)
    tacs_aim.add_variable(variable=thick_DV)

load = caps2tacs.GridForce(name="load1", caps_load="plate", direction=[0,0,-1.], magnitude=1.0E3)    
tacs_aim.add_load(load=load)

tacs_aim.setup_aim()
egads_aim.set_mesh()


pytacs_function = caps2tacs.MassStressTransient(t0=0.0, tf=10.0, num_steps=3, amplitude=lambda t : np.sin(0.5*t))

# start a caps tacs main problem
caps_tacs = caps2tacs.CapsTacs(name="transient_panel", tacs_aim=tacs_aim, egads_aim=egads_aim, pytacs_function=pytacs_function, compute_gradients=False)

caps_tacs.analysis()

# convert all f5 to vtk in analysis dir 
for filename in os.listdir(tacs_aim.analysis_dir):
    filepath = os.path.join(tacs_aim.analysis_dir, filename)
    os.system(f"~/git/tacs/extern/f5tovtk/f5tovtk {filepath}")

print(f"paraview ./capsStruct/Scratch/tacs/mass_stress_0_000..vtk")





