import capsManager

panel_problem = capsManager.CapsFluid.default(csmFile="diamondWing.csm")
pointwise_aim = panel_problem.pointwiseAim
fun3d_aim = panel_problem.fun3dAim

pointwise_aim.set_mesh()
pointwise_aim.run_pointwise()