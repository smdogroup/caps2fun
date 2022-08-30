
import os
from capsManager.fun3d.flow_settings import FlowSettings, MotionSettings
from funtofemManager.analysis import Analysis
__all__ = ["Fun3dAim"]


from typing import TYPE_CHECKING
import pyCAPS
from capsManager.fun3d.flow_settings import FlowSettings


class Fun3dAim:
    def __init__(self, caps_problem:pyCAPS.Problem, flow_settings:FlowSettings, motion_settings:MotionSettings):
        self._aim = caps_problem.analysis.create(aim = "fun3dAIM",
                                    name = "fun3d")
        self._flow_settings = flow_settings
        self._motion_settings = motion_settings

        self._is_setup = True

    @property
    def is_setup(self) -> bool:
        return self._is_setup
                                    
    @property
    def aim(self):
        return self._aim

    @property
    def flow_settings(self) -> FlowSettings:
        return self._flow_settings

    @property
    def motion_settings(self) -> MotionSettings:
        return self._motion_settings

    @property
    def analysis_type(self) -> str:
        return self._analysis_type

    @property
    def analysis_dir(self) -> str:
        return self.aim.analysisDir  

    @property
    def flow_directory(self) -> str:
        return os.path.join(self.analysis_dir, "Flow")

    @property
    def adjoint_directory(self) -> str:
        return os.path.join(self.analysis_dir, "Adjoint")    

    @property
    def project_name(self) -> str:
        return self.aim.projName

    def pre_analysis(self):
        self.aim.preAnalysis()

    def post_analysis(self):
        self.aim.postAnalysis()

    
        