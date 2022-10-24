
__all__ = ["Fun3dAim"]


from typing import TYPE_CHECKING
import pyCAPS, os
from capsManager.fun3d.fun3d_settings import FlowSettings, FluidMeshSettings, MotionSettings, FluidSolverSettings


class Fun3dAim:
    def __init__(self, 
        caps_problem:pyCAPS.Problem, 
        flow_settings:FlowSettings=None, 
        motion_settings:MotionSettings=None, 
        fluid_mesh_settings:FluidMeshSettings=None,
        fluid_solver_settings:FluidSolverSettings=None,
        build_complex:bool=False,
        write_vtk:bool=False
    ):
        self._aim = caps_problem.analysis.create(aim = "fun3dAIM",
                                    name = "fun3d")
        self._flow_settings = flow_settings
        self._motion_settings = motion_settings
        self._fluid_mesh_settings = fluid_mesh_settings if fluid_mesh_settings is not None else FluidMeshSettings()
        self._fluid_solver_settings = fluid_solver_settings if fluid_solver_settings is not None else FluidSolverSettings()
        self._build_complex = build_complex

        # set to not overwrite fun3d nml analysis
        self.aim.input.Overwrite_NML = False
        #fun3d design sensitivities settings
        self.aim.input.Design_SensFile = True
        self.aim.input.Design_Sensitivity = True

        self._write_vtk = write_vtk

    @property
    def is_setup(self) -> bool:
        return self._flow_settings is not None and self._motion_settings is not None
                                    
    @property
    def aim(self):
        return self._aim

    @property
    def flow_settings(self) -> FlowSettings:
        return self._flow_settings

    @flow_settings.setter
    def flow_settings(self, new_settings:FlowSettings):
        self._flow_settings = new_settings

    @property
    def motion_settings(self) -> MotionSettings:
        return self._motion_settings

    @motion_settings.setter
    def motion_settings(self, new_settings:MotionSettings):
        self._motion_settings = new_settings

    @property
    def fluid_mesh_settings(self) -> FluidMeshSettings:
        return self._fluid_mesh_settings

    @fluid_mesh_settings.setter
    def fluid_mesh_settings(self, new_settings:FluidMeshSettings):
        self._fluid_mesh_settings = new_settings

    @property
    def fluid_solver_settings(self) -> FluidMeshSettings:
        return self._fluid_solver_settings

    @fluid_solver_settings.setter
    def fluid_solver_settings(self, new_settings:FluidSolverSettings):
        self._fluid_solver_settings = new_settings

    @property
    def build_complex(self) -> bool:
        return self._build_complex

    @build_complex.setter
    def build_complex(self, new_bool:bool):
        self._build_complex = new_bool
        if new_bool is True:
            self._print_perturb = True

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
        return self.aim.input.Proj_Name

    def pre_analysis(self):
        self.aim.preAnalysis()

    def post_analysis(self):
        self.aim.postAnalysis()

    def set_boundary_condition(self, boundary_condition_dict:dict):
        self.aim.input.Boundary_Condition = boundary_condition_dict

    def write(self):
        """
        generate the writer classes for fun3d.nml and moving_body.input files and then write them
        """
        from capsManager.fun3d.namelist_writer import Fun3dNamelistWriter, MovingBodyInputWriter, print_perturb_input
        Fun3dNamelistWriter(fun3d_aim=self).write()
        MovingBodyInputWriter(fun3d_aim=self).write()
        if self.build_complex:
            print_perturb_input(path=self.flow_directory)
        