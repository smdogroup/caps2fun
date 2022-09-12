import caps2fun
# if not using run.pbs, load pointwise first
# module load pointwise/18.5R1

caps_fluid = caps2fun.CapsFluid.default(csmFile="wedge_fluid.csm")
pointwise_aim = caps_fluid.pointwiseAim
fun3d_aim = caps_fluid.fun3dAim

pointwise_aim.set_mesh(
        wall_spacing=0.01, 
        domain_max_layers=15,
        domain_growth_rate=1.75,
        block_boundary_decay=0.5,
        block_max_skew_angle=170.0,
        block_full_layers=1,
        block_max_layers=100,
        inviscid=True
)
fun3d_aim.flow_settings = caps2fun.FlowSettings(
        flow_type="inviscid",
        mach_number=0.3,
        angle_of_attack=0.0,
        reynolds_number=1e4,
        temperature_ref=300.0,
        ref_area=1.0,
        num_steps=20,
        freeze_limiter_iteration=None
        )
fun3d_aim.motion_settings = caps2fun.MotionSettings(body_name="wedge")
fun3d_aim.build_complex = True

caps_fun3d = caps2fun.CapsFun3d(pointwise_aim=pointwise_aim, fun3d_aim=fun3d_aim)
caps_fun3d.build_mesh()
caps_fun3d.prepare_fun3d()