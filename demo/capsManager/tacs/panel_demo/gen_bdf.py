from curses import panel
import capsWrapper

panel_problem = capsWrapper.CapsStruct.default(csmFile="panel.csm")
tacs_aim = panel_problem.tacsAim
egads_aim = panel_problem.egadsAim

madeupium = capsWrapper.Isotropic.madeupium()
tacs_aim.add_material(material=madeupium)


constraint = capsWrapper.ZeroConstraint(name="fixEdge", caps_constraint="edge")
tacs_aim.add_constraint(constraint=constraint)

load = capsWrapper.GridForce(name="load1", caps_load="plate", direction=[0,0,-1.], magnitude=1.0E3)    
tacs_aim.add_load(load=load)

thick_DV = capsWrapper.ThicknessVariable(name="thick", caps_group="plate", value=0.01, material=madeupium)
tacs_aim.add_variable(variable=thick_DV)

despmtrs = ["plateLength", "plateWidth"]
for despmtr in despmtrs:
    shape_var = capsWrapper.ShapeVariable(name=despmtr)
    tacs_aim.add_variable(variable=shape_var)

tacs_aim.setup_aim()

built_tacs_aim = tacs_aim.aim
built_egads_aim = egads_aim.aim

built_tacs_aim.input["Mesh"].link(built_egads_aim.output["Surface_Mesh"])

# run the tacs aim preanalysis
built_tacs_aim.preAnalysis()



