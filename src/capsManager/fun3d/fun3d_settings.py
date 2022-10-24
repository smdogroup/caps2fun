

__all__ = ["FlowSettings", "MotionSettings", "FluidMeshSettings", "FluidSolverSettings"]

from typing import TYPE_CHECKING, List

class FlowSettings:
    """
    Fun3d flow settings class
    """
    FLOW_TYPES = ["inviscid", "laminar", "turbulent"]
    def __init__(self, 
        flow_type:str="inviscid",
        mach_number:float=0.3,
        angle_of_attack:float=0.0,
        reynolds_number:float=1e6,
        temperature_ref:float=300.0,
        ref_area:float=1.0,
        num_steps:int=10.0,
        freeze_limiter_iteration:int=None
        ):
        assert(flow_type in FlowSettings.FLOW_TYPES)
        self._flow_type = flow_type
        self._mach_number = mach_number
        self._angle_of_attack = angle_of_attack
        self._reynolds_number = reynolds_number
        self._temperature_ref = temperature_ref
        self._ref_area = ref_area
        self._num_steps = num_steps
        if freeze_limiter_iteration is None or freeze_limiter_iteration > num_steps:
            freeze_limiter_iteration = int(5/6 * num_steps)
        self._freeze_limiter_iteration = freeze_limiter_iteration

    @property
    def mach_number(self) -> float:
        return self._mach_number

    @property
    def angle_of_attack(self) -> float:
        return self._angle_of_attack

    @property
    def flow_type(self) -> str:
        return self._flow_type

    @property
    def reference_temperature(self) -> float:
        return self._temperature_ref

    @property
    def reynolds_number(self) -> float:
        return self._reynolds_number

    @property
    def reference_area(self) -> float:
        return self._ref_area

    @property
    def num_steps(self) -> int:
        return self._num_steps

    @property
    def freeze_limiter_iteration(self) -> int:
        return self._freeze_limiter_iteration

class FluidSolverSettings:
    def __init__(self,
        schedule_iteration:List[int]=[1,50],
        schedule_cfl:List[float]=[200,200],
        time_accuracy:str='steady',
        sub_iterations:int=0,
        schedule_cfl_turb:List[float] = [50.0,50.0]
    ):
        self.schedule_iteration = schedule_iteration
        self.schedule_cfl = schedule_cfl
        self.time_accuracy = time_accuracy
        self.sub_iterations = sub_iterations
        self.schedule_cfl_turb = schedule_cfl_turb

class FluidMeshSettings:
    def __init__(self,
        num_search:int=200,
        tolerance:float=1.e-14,
        substeps:int=1,
        from_initial:bool=True,
        use_substeps:bool=False,
        elasticity_const:int=1,
        elasticity_exponent:float=1.0,
        num_restarts:int=1,
        poisson_ratio:float=0.0
    ):
        self.num_search = num_search
        self.tolerance = tolerance
        self.substeps = substeps
        self.from_initial = from_initial
        self.use_substeps = use_substeps
        self.elasticity_const = int(elasticity_const)
        self.elasticity_exponent = elasticity_exponent
        self.num_restarts = num_restarts
        self.poisson_ratio = poisson_ratio

class MotionSettings:
    def __init__(self,
        body_name:str, 
        num_bodies:int=1,
        motion_style:str="deform",
        defining_boundary:int=1
    ):
        """
        TODO : currently only supports one body, one boundary, update later
        """
        assert(motion_style in ["deform", "rigid", "deform+rigid", "rigid+deform"])
        self.body_name = body_name
        self.num_bodies = num_bodies
        self.defining_boundary = defining_boundary
        self.motion_style = motion_style