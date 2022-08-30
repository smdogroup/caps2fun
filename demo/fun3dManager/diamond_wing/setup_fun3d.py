import caps2fun
from capsManager import fun3d

caps_fluid = caps2fun.CapsFluid.default(csmFile="diamondWing.csm")
pointwise_aim = caps2fun.PointwiseAim(caps_problem=caps_fluid)

flow_settings = caps2fun.FlowSettings()
motion_settings = caps2fun.MotionSettings(body_name="diamondWing")
fun3d_aim = caps2fun.Fun3dAim(caps_problem=caps_fluid, flow_settings=flow_settings, motion_settings=motion_settings)

caps_fun3d = caps2fun.CapsFun3d(pointwise_aim=pointwise_aim, fun3d_aim=fun3d_aim)
caps_fun3d.build_mesh()
caps_fun3d.prepare_fun3d()