

from typing import TYPE_CHECKING, List
import f90nml, os
from capsManager.fun3d.fun3d_aim import Fun3dAim

class Fun3dNamelistWriter:
    def __init__(self, fun3d_aim:Fun3dAim):

        self._fun3d_aim = fun3d_aim
        self._fun3d_nml = f90nml.Namelist()

        self._set_namelist()

    @property
    def aim(self):
        return self._fun3d_aim

    @property
    def fun3d_nml(self):
        return self._fun3d_nml

    def _set_namelist(self):
        
        #project section
        self.fun3d_nml["project"] = f90nml.Namelist()
        self.fun3d_nml["project"]["project_rootname"] = self.aim.project_name
        #self.fun3dAim.input.Proj_Name = self.config["mesh_style"]

        #governing equation section
        self.fun3d_nml["governing_equations"] = f90nml.Namelist()

        #apply fun3d analysis type
        #options = "inviscid", "laminar", "turbulent"
        if (self.aim.flow_settings.flow_type == "inviscid"):
            self.fun3d_nml["governing_equations"]["viscous_terms"] = "inviscid"
        elif (self.aim.flow_settings.flow_type == "laminar"):
            self.fun3d_nml["governing_equations"]["eqn_type"] = "compressible"
            self.fun3d_nml["governing_equations"]["viscous_terms"] = "laminar"
        elif (self.aim.flow_settings.flow_type == "turbulent"):
            raise AssertionError("Turbulent not supported yet")

        #raw grid section
        self.fun3d_nml["raw_grid"] = f90nml.Namelist()
        #self.fun3dAim.input.Mesh_ASCII_Flag = False
        self.fun3d_nml["raw_grid"]["grid_format"] = "aflr3"
        self.fun3d_nml["raw_grid"]["data_format"] = "default"
        self.fun3d_nml["raw_grid"]["swap_yz_axes"] = False

        #reference physical properties section
        self.fun3d_nml["reference_physical_properties"] = f90nml.Namelist()
        self.fun3d_nml["reference_physical_properties"]["mach_number"] = self.aim.flow_settings.mach_number
        self.fun3d_nml["reference_physical_properties"]["angle_of_attack"] = self.aim.flow_settings.angle_of_attack
        self.fun3d_nml["reference_physical_properties"]["reynolds_number"] = self.aim.flow_settings.reynolds_number
        self.fun3d_nml["reference_physical_properties"]["temperature"] = self.aim.flow_settings.reference_temperature
        self.fun3d_nml["reference_physical_properties"]["temperature_units"] = "Kelvin"
        
        #inviscid flux method section
        self.fun3d_nml["inviscid_flux_method"] = f90nml.Namelist()
        self.fun3d_nml["inviscid_flux_method"]["flux_construction"] = "roe"
        self.fun3d_nml["inviscid_flux_method"]["flux_limiter"] = "hminmod"
        self.fun3d_nml["inviscid_flux_method"]["smooth_limiter_coeff"] = 1.0
        self.fun3d_nml["inviscid_flux_method"]["freeze_limiter_iteration"] = self.aim.flow_settings.freeze_limiter_iteration

        #nonlinear solver parameters section
        self.fun3d_nml["nonlinear_solver_parameters"] = f90nml.Namelist()
        self.fun3d_nml["nonlinear_solver_parameters"]["schedule_iteration"] = [1, 80]
        self.fun3d_nml["nonlinear_solver_parameters"]["schedule_cfl"] = [2, 100]
        if (self.aim.flow_settings.flow_type == "inviscid"):
            self.fun3d_nml["nonlinear_solver_parameters"]["time_accuracy"] = "steady"
            self.fun3d_nml["nonlinear_solver_parameters"]["time_step_nondim"] = 0.1
            self.fun3d_nml["nonlinear_solver_parameters"]["subiterations"] = 0
            self.fun3d_nml["nonlinear_solver_parameters"]["schedule_cflturb"] = [50.0,50.0]


        #force moment integ properties section
        self.fun3d_nml["force_moment_integ_properties"] = f90nml.Namelist()
        self.fun3d_nml["force_moment_integ_properties"]["area_reference"] = self.aim.flow_settings.reference_area

        #code run and control section
        self.fun3d_nml["code_run_control"] = f90nml.Namelist()
        self.fun3d_nml["code_run_control"]["steps"] = self.aim.flow_settings.num_steps
        self.fun3d_nml["code_run_control"]["stopping_tolerance"] = 1.0e-15
        self.fun3d_nml["code_run_control"]["restart_write_freq"] = 1000
        self.fun3d_nml["code_run_control"]["restart_read"] = "off"

        #global flow_settings
        self.fun3d_nml["global"] = f90nml.Namelist()
        self.fun3d_nml["global"]["moving_grid"] = True
        self.fun3d_nml["global"]["volume_animation_freq"] = -1
        self.fun3d_nml["global"]["boundary_animation_freq"] = -1

        #mesh elasticity flow_settings
        self.fun3d_nml["elasticity_gmres"] = f90nml.Namelist()
        self.fun3d_nml["elasticity_gmres"]["algebraic_mesh_deform"] = False #default False
        self.fun3d_nml["elasticity_gmres"]["nsearch"] = 200 #default 50
        self.fun3d_nml["elasticity_gmres"]["tol"] = 1.e-14
        self.fun3d_nml["elasticity_gmres"]["deformation_substeps"] = 1 #default 1, other value 5
        self.fun3d_nml["elasticity_gmres"]["deform_from_initial_mesh"] = True #default true
        self.fun3d_nml["elasticity_gmres"]["use_substeps_each_step"] = False #default False, need to turn to True if deformation_substeps used
        self.fun3d_nml["elasticity_gmres"]["elasticity"] = 1 #default 1, option 2
        self.fun3d_nml["elasticity_gmres"]["elasticity_exponent"] = 1.0 #default 1.0, change to 2.0 if needed
        self.fun3d_nml["elasticity_gmres"]["nrestarts"] = 1
        self.fun3d_nml["elasticity_gmres"]["poisson_ratio"] = 0.0 #default 0.0
        
        #massoud output flow_settings
        self.fun3d_nml["massoud_output"] = f90nml.Namelist()
        self.fun3d_nml["massoud_output"]["funtofem_include_skin_friction"] = False

        #volume output variables
        if self.aim._write_vtk:
            self.fun3d_nml["volume_output_variables"] = f90nml.Namelist()
            self.fun3d_nml["volume_output_variables"]["export_to"] = "vtk"
            self.fun3d_nml["volume_output_variables"]["x"] = False
            self.fun3d_nml["volume_output_variables"]["y"] = False
            self.fun3d_nml["volume_output_variables"]["z"] = False
            self.fun3d_nml["volume_output_variables"]["temperature"] = True
            self.fun3d_nml["volume_output_variables"]["mach"] = True
            self.fun3d_nml["volume_output_variables"]["p"] = True

        #boundary output variables
        if self.aim._write_vtk:
            self.fun3d_nml["boundary_output_variables"] = f90nml.Namelist()
            #boundary list indexes probably auto set from fun3dAim
            self.fun3d_nml["boundary_output_variables"]["number_of_boundaries"] = -1
            self.fun3d_nml["boundary_output_variables"]["boundary_list"] = "1-2"
            self.fun3d_nml["boundary_output_variables"]["temperature"] = True
            self.fun3d_nml["boundary_output_variables"]["mach"] = True
            self.fun3d_nml["boundary_output_variables"]["p"] = True

    def write(self):
        """
        writes the fun3d namelist file 'fun3d.nml' in the fun3d Flow directory
        """
        fun3d_nml_file = os.path.join(self.aim.flow_directory, "fun3d.nml")
        self.fun3d_nml.write(fun3d_nml_file, force=True)

class MovingBodyInputWriter:
    """
    Fun3d class to write moving_body.input file in fun3d flow directory
    """
    def __init__(self, fun3d_aim:Fun3dAim):
        self._fun3d_aim = fun3d_aim
        self._moving_body_input = f90nml.Namelist()

        self._set_namelist()

    @property
    def aim(self) -> Fun3dAim:
        return self._fun3d_aim

    def _set_namelist(self):
        #moving body flow_settings for funtofem to fun3d
        bndryArray = self.aim.motion_settings.boundary_indices
        bndryArray = list(bndryArray)

        #body definitions
        self._moving_body_input["body_definitions"] = f90nml.Namelist()
        self._moving_body_input["body_definitions"]["n_moving_bodies"] = self.aim.motion_settings.num_bodies
        self._moving_body_input["body_definitions"]["body_name"] = [self.aim.motion_settings.body_name]
        self._moving_body_input["body_definitions"]["parent_name"] = [""] # '' means motion relative to inertial ref frame
        self._moving_body_input["body_definitions"]["n_defining_bndry"] = [self.aim.motion_settings.num_boundaries] #number of boundaries that define this body
        self._moving_body_input["body_definitions"]["defining_bndry(1,1)"] = 2 #index 1: boundary number index 2: body number
        self._moving_body_input["body_definitions"]["motion_driver"] = ["funtofem"] #tells fun3d to use motion inputs from python
        self._moving_body_input["body_definitions"]["mesh_movement"] = [self.aim.motion_settings.motion_style] #can use 'rigid', 'deform', 'rigid+deform' with funtofem interface

        self._setup = True

    def write(self):
        moving_body_input_file = os.path.join(self.aim.flow_directory, "moving_body.input")
        self._moving_body_input.write(moving_body_input_file, force=True)

perturb_file = """   PERTURB   EPSILON GRIDPOINT
        0    1e-30       666

0 = No perturbation
1 = Mach number
2 = Alpha
3 = Shape
4 = x-rotation rate
5 = y-rotation rate
6 = z-rotation rate
7 = Grid point x
8 = Grid point y
9 = Grid point z
10 = Yaw
11 = error transport (truncation error)
12 = RCS jet plenum pressure, p0
100+ = add an imaginary source term to equation
        PERTURB-100 of node GRIDPOINT
        (to verify the adjoint lambda value)
"""
def print_perturb_input(path:str, comm=None):
    if comm is None or comm.Get_rank() == 0:
        filepath = os.path.join(path, "perturb.input")
        file_hdl = open(filepath, "w")
        file_hdl.write(perturb_file)
        file_hdl.close()
