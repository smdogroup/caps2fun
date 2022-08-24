import caps2tacs
from pyoptsparse import SLSQP, Optimization

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
egads_aim.set_mesh(edge_pt_min=20, edge_pt_max=40, global_mesh_size=0.01, max_surf_offset=0.01, max_dihedral_angle=5)

# make a pytacs function
pytacs_function = caps2tacs.MassStress()

caps_tacs = caps2tacs.CapsTacs(
    name="naca_wing_struct", tacs_aim=tacs_aim, 
    egads_aim=egads_aim, pytacs_function=pytacs_function, view_plots=True
    )


def initial_stress(relative:bool=True):
    if relative:
        return caps_tacs.function(function_name="ks_vmfailure")
    else:
        return 1.0

names = ["mass", "stress"]

def objfunc(xdict):
    funcs = {}
    fail = False
    try:
        funcs["mass"] = caps_tacs.function("mass", xdict)
        funcs["stress"] = caps_tacs.function("ks_vmfailure", xdict)
    except:
        fail = True
        for name in names:
            funcs[name] = None

    print(funcs)
    #sys.exit()

    return funcs, fail

def objsens(xdict, funcs):
    grads = {}
    fail = False
    try:
        grads["mass"] = caps_tacs.gradient("mass", xdict)
        grads["stress"] = caps_tacs.gradient("ks_vmfailure", xdict)
    except:
        for name in names:
            grads[name] = None
        fail = True

    print(grads)

    return grads, fail

opt_problem = Optimization(name="purdue_p_structure", objFun=objfunc)
opt_problem.addObj(name="mass", scale=0.001)
opt_problem.addCon(name="stress", upper=initial_stress(relative=False), scale=1.0)

for dv in tacs_aim.thickness_design_variables:
    opt_problem.addVar(name=dv.name, lower=0.0001, upper=0.50, value=dv.value, scale=100.0)

for shape_var in tacs_aim.shape_design_variables:
    print(shape_var.name, shape_var.value)
    opt_problem.addVar(name=shape_var.name, value=shape_var.value, lower=0.5*shape_var.value, upper=2.0*shape_var.value, scale=100.0)

slsqp_method = SLSQP(options={
    "IPRINT" : -1,
    "MAXIT" : 10
})

solution = slsqp_method(opt_problem, sens=objsens)

# rst begin check
# Check Solution
print(solution)




