import caps2tacs

panel_problem = caps2tacs.CapsStruct.default(csmFile="simple_diamond.csm")
tacs_aim = panel_problem.tacsAim
egads_aim = panel_problem.egadsAim

aluminum = caps2tacs.Isotropic.aluminum()
tacs_aim.add_material(material=aluminum)

# add wing root constraint
constraint = caps2tacs.ZeroConstraint(name="fixRoot", caps_constraint="root")
tacs_aim.add_constraint(constraint=constraint)

Tconstraint = caps2tacs.TemperatureConstraint(name="tempRoot", caps_constraint="root", temperature=300)
tacs_aim.add_constraint(constraint=Tconstraint)

thick_DV = caps2tacs.ThicknessVariable(name="thick",caps_group="shell",value=0.2,material=aluminum)
tacs_aim.add_variable(thick_DV)

tacs_aim.setup_aim(auto_shape_variables=True)
egads_aim.set_mesh(edge_pt_min=20, edge_pt_max=30, global_mesh_size=0.10, max_surf_offset=0.01, max_dihedral_angle=5)

caps_tacs = caps2tacs.CapsTacs(
    name="simple_diamond", tacs_aim=tacs_aim, 
    egads_aim=egads_aim, pytacs_function=None,
    )

# build the mesh
caps_tacs.tacs_aim.pre_analysis()