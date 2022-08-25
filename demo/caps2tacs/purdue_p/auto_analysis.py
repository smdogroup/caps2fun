import caps2tacs

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

load = caps2tacs.GridForce(name="load1", caps_load="front", direction=[0,0,1.0], magnitude=1.0E3)    
tacs_aim.add_load(load=load)

tacs_aim.setup_aim()
egads_aim.set_mesh(edge_pt_min=20, edge_pt_max=30, global_mesh_size=0.05)

# make a pytacs function
pytacs_function = caps2tacs.MassStress()

# start a caps tacs main problem
caps_tacs = caps2tacs.CapsTacs(tacs_aim=tacs_aim, egads_aim=egads_aim, pytacs_function=pytacs_function)

mass_grad,mass_grad_dict = caps_tacs.gradient(function_name="mass")
print(f"mass gradient dict = {mass_grad_dict}\n")

stress_grad,stress_grad_dict = caps_tacs.gradient(function_name="ks_vmfailure")
print(f"stress gradient dict = {stress_grad_dict}\n")



