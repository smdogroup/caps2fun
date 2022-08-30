import caps2tacs
import os
import numpy as np

# all the capsWrapper package modules/classes are also available 
# in caps2tacs since caps2tacs package imports the whole capsWrapper package\

panel_problem = caps2tacs.CapsStruct.default(csmFile="plane.csm")
tacs_aim = panel_problem.tacsAim
egads_aim = panel_problem.egadsAim

aluminum = caps2tacs.Isotropic.aluminum()
tacs_aim.add_material(material=aluminum)

constraint = caps2tacs.ZeroConstraint(name="fixEdge", caps_constraint="rear")
tacs_aim.add_constraint(constraint=constraint)

load = caps2tacs.GridForce(name="load1", caps_load="wing", direction=[0,-1.0,0.0], magnitude=1.0E2)    
tacs_aim.add_load(load=load)

thick_idx = 0

# caps_groups = [
#     wing:OML, 
#     fuselage:OML, 
#     wing:spar3, 
#     wing:rib2, 
#     fuselage:station3, 
#     fuselage:station4, 
#     wing:spar2, 
#     wing:rib3, 
#     fuselage:station2, 
#     wing:spar1, 
#     fuselage:station5, 
#     wing:rib4, 
#     nose:OML, 
#     nose:station1, fuselage:station6, wing:rib5, nose:station2, fuselage:station7, wing:rib6, nose:station3, fuselage:station8, wing:rib7, fuselage:station9, htail:OML, wing:rib8, fuselage:station10, htail:spar1, htail:rib1, wing:rib9, htail:spar2, htail:rib2, wing:rib10, htail:rib3, wing:rib11, htail:rib4, wing:rib12, htail:rib5, wing:rib13, htail:rib6, wing:rib14, htail:rib7, wing:rib15, htail:rib8, wing:rib16,
# ]

section = "wing"
caps_group_suffixes = [f"rib{idx+1}" for idx in range(16)] + [f"spar{idx+1}" for idx in range(3)] + ["OML"]
for caps_group_suffix in caps_group_suffixes:
    if not(caps_group_suffix in ["rib1"]):
        if "rib" in caps_group_suffix:
            thickness = 0.025
        elif "spar" in caps_group_suffix:
            thickness = 0.03
        else:
            thickness = 0.01
        thick_DV = caps2tacs.ThicknessVariable(
            name=f"thick{thick_idx}", caps_group=f"{section}:{caps_group_suffix}", 
            value=thickness, material=aluminum
            )
        tacs_aim.add_variable(variable=thick_DV)
        thick_idx += 1

section = "fuselage"
caps_group_suffixes = [f"station{idx+1}" for idx in range(10)] + ["OML"]
for caps_group_suffix in caps_group_suffixes:
    if not(caps_group_suffix in ["station1"]):
        if "station" in caps_group_suffix:
            thickness = 0.03
        else:
            thickness = 0.01
        thick_DV = caps2tacs.ThicknessVariable(
            name=f"thick{thick_idx}", caps_group=f"{section}:{caps_group_suffix}", 
            value=thickness, material=aluminum
            )
        tacs_aim.add_variable(variable=thick_DV)
        thick_idx += 1

section = "nose"
caps_group_suffixes = [f"station{idx+1}" for idx in range(3)] + ["OML"]
for caps_group_suffix in caps_group_suffixes:
    if not(caps_group_suffix in []):
        if "station" in caps_group_suffix:
            thickness = 0.03
        else:
            thickness = 0.01
        thick_DV = caps2tacs.ThicknessVariable(
            name=f"thick{thick_idx}", caps_group=f"{section}:{caps_group_suffix}", 
            value=thickness, material=aluminum
            )
        tacs_aim.add_variable(variable=thick_DV)
        thick_idx += 1

section = "htail"
caps_group_suffixes = [f"rib{idx+1}" for idx in range(8)] + [f"spar{idx+1}" for idx in range(2)] + ["OML"]
for caps_group_suffix in caps_group_suffixes:
    if not(caps_group_suffix in []):
        if "rib" in caps_group_suffix:
            thickness = 0.025
        elif "spar" in caps_group_suffix:
            thickness = 0.03
        else:
            thickness = 0.01
        thick_DV = caps2tacs.ThicknessVariable(
            name=f"thick{thick_idx}", caps_group=f"{section}:{caps_group_suffix}", 
            value=thickness, material=aluminum
            )
        tacs_aim.add_variable(variable=thick_DV)
        thick_idx += 1


tacs_aim.setup_aim()
egads_aim.set_mesh()

pytacs_function = caps2tacs.MassStress()

# start a caps tacs main problem
caps_tacs = caps2tacs.CapsTacs(name="static_aircraft", tacs_aim=tacs_aim, egads_aim=egads_aim, pytacs_function=pytacs_function, 
compute_gradients=False, view_plots=True)

caps_tacs.analysis()

