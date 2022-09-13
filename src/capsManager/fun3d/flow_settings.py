

__all__ = ["FlowSettings", "MotionSettings"]

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

class MotionSettings:
    def __init__(self,
        body_name:str, 
        num_bodies:int=1,
        num_boundaries:int=1,
        boundary_indices:List[int]=[2],
        motion_style:str="deform"
    ):
        assert(motion_style in ["deform", "rigid", "deform+rigid", "rigid+deform"])
        self._body_name = body_name
        self._num_bodies = num_bodies
        self._num_boundaries = num_boundaries
        self._boundary_indices = boundary_indices
        self._motion_style = motion_style

    @property
    def body_name(self):
        return self._body_name

    @property
    def num_bodies(self):
        return self._num_bodies

    @property
    def num_boundaries(self):
        return self._num_boundaries

    @property
    def boundary_indices(self):
        return self._boundary_indices

    @property
    def motion_style(self):
        return self._motion_style