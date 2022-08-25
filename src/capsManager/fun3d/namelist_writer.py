

class Fun3dNamelistWriter:
    def __init__(self):
        """
        # namelist and mapbc settings
        self.capsFluid.analysis["fun3d"].input.Overwrite_NML = False
        self.fun3dnml = f90nml.Namelist()
        
        #project section
        self.fun3dnml["project"] = f90nml.Namelist()
        self.fun3dnml["project"]["project_rootname"] = self.fun3dAim.input.Proj_Name
        #self.fun3dAim.input.Proj_Name = self.config["mesh_style"]

        #governing equation section
        self.fun3dnml["governing_equations"] = f90nml.Namelist()

        #apply fun3d analysis type
        #options = "inviscid", "laminar", "turbulent"
        if (self.config["fun3d_analysis_type"] == "inviscid"):
            self.fun3dnml["governing_equations"]["viscous_terms"] = "inviscid"
        elif (self.config["fun3d_analysis_type"] == "laminar"):
            self.fun3dnml["governing_equations"]["eqn_type"] = "compressible"
            self.fun3dnml["governing_equations"]["viscous_terms"] = "laminar"
        elif (self.config["fun3d_analysis_type"] == "turbulent"):
            pass

        #raw grid section
        self.fun3dnml["raw_grid"] = f90nml.Namelist()
        #self.fun3dAim.input.Mesh_ASCII_Flag = False
        self.fun3dnml["raw_grid"]["grid_format"] = "aflr3"
        self.fun3dnml["raw_grid"]["data_format"] = "default"
        self.fun3dnml["raw_grid"]["swap_yz_axes"] = False

        #reference physical properties section
        self.fun3dnml["reference_physical_properties"] = f90nml.Namelist()
        self.fun3dnml["reference_physical_properties"]["mach_number"] = self.config["mach"]
        self.fun3dnml["reference_physical_properties"]["angle_of_attack"] = self.config["AOA"]
        self.fun3dnml["reference_physical_properties"]["reynolds_number"] = self.config["Re"]
        self.fun3dnml["reference_physical_properties"]["temperature"] = self.config["temperature"]
        self.fun3dnml["reference_physical_properties"]["temperature_units"] = "Kelvin"
        #self.capsFluid.analysis["fun3d"].input.Alpha = 1.0
        #self.capsFluid.analysis["fun3d"].input.Mach = 0.5
        #self.capsFluid.analysis["fun3d"].input.Re = 35e6
        
        #inviscid flux method section
        self.fun3dnml["inviscid_flux_method"] = f90nml.Namelist()
        self.fun3dnml["inviscid_flux_method"]["flux_construction"] = "roe"
        self.fun3dnml["inviscid_flux_method"]["flux_limiter"] = "hminmod"
        self.fun3dnml["inviscid_flux_method"]["smooth_limiter_coeff"] = 1.0
        self.fun3dnml["inviscid_flux_method"]["freeze_limiter_iteration"] = int(5.0/6 * self.config["nsteps"])

        #nonlinear solver parameters section
        self.fun3dnml["nonlinear_solver_parameters"] = f90nml.Namelist()
        self.fun3dnml["nonlinear_solver_parameters"]["schedule_iteration"] = [1, 80]
        self.fun3dnml["nonlinear_solver_parameters"]["schedule_cfl"] = [2, 100]
        if (not(self.config["fun3d_analysis_type"] == "inviscid")):
            self.fun3dnml["nonlinear_solver_parameters"]["time_accuracy"] = "steady"
            self.fun3dnml["nonlinear_solver_parameters"]["time_step_nondim"] = 0.1
            self.fun3dnml["nonlinear_solver_parameters"]["subiterations"] = 0
            self.fun3dnml["nonlinear_solver_parameters"]["schedule_cflturb"] = [50.0,50.0]


        #force moment integ properties section
        self.fun3dnml["force_moment_integ_properties"] = f90nml.Namelist()
        self.fun3dnml["force_moment_integ_properties"]["area_reference"] = self.capsStruct.geometry.despmtr["area"].value / 2

        #code run and control section
        self.fun3dnml["code_run_control"] = f90nml.Namelist()
        self.fun3dnml["code_run_control"]["steps"] = self.config["nsteps"]
        self.fun3dnml["code_run_control"]["stopping_tolerance"] = 1.0e-15
        self.fun3dnml["code_run_control"]["restart_write_freq"] = 1000
        self.fun3dnml["code_run_control"]["restart_read"] = "off"

        #global settings
        self.fun3dnml["global"] = f90nml.Namelist()
        self.fun3dnml["global"]["moving_grid"] = True
        self.fun3dnml["global"]["volume_animation_freq"] = -1
        self.fun3dnml["global"]["boundary_animation_freq"] = -1

        #mesh elasticity settings
        self.fun3dnml["elasticity_gmres"] = f90nml.Namelist()
        self.fun3dnml["elasticity_gmres"]["algebraic_mesh_deform"] = False #default False
        self.fun3dnml["elasticity_gmres"]["nsearch"] = 200 #default 50
        self.fun3dnml["elasticity_gmres"]["tol"] = 1.e-14
        self.fun3dnml["elasticity_gmres"]["deformation_substeps"] = 1 #default 1, other value 5
        self.fun3dnml["elasticity_gmres"]["deform_from_initial_mesh"] = True #default true
        self.fun3dnml["elasticity_gmres"]["use_substeps_each_step"] = False #default False, need to turn to True if deformation_substeps used
        self.fun3dnml["elasticity_gmres"]["elasticity"] = 1 #default 1, option 2
        self.fun3dnml["elasticity_gmres"]["elasticity_exponent"] = 1.0 #default 1.0, change to 2.0 if needed
        self.fun3dnml["elasticity_gmres"]["nrestarts"] = 1
        self.fun3dnml["elasticity_gmres"]["poisson_ratio"] = 0.0 #default 0.0
        
        #massoud output settings
        self.fun3dnml["massoud_output"] = f90nml.Namelist()
        self.fun3dnml["massoud_output"]["funtofem_include_skin_friction"] = False

        #volume output variables
        self.fun3dnml["volume_output_variables"] = f90nml.Namelist()
        self.fun3dnml["volume_output_variables"]["export_to"] = "vtk"
        self.fun3dnml["volume_output_variables"]["x"] = False
        self.fun3dnml["volume_output_variables"]["y"] = False
        self.fun3dnml["volume_output_variables"]["z"] = False
        self.fun3dnml["volume_output_variables"]["temperature"] = True
        self.fun3dnml["volume_output_variables"]["mach"] = True
        self.fun3dnml["volume_output_variables"]["p"] = True

        #boundary output variables
        self.fun3dnml["boundary_output_variables"] = f90nml.Namelist()
        #boundary list indexes probably auto set from fun3dAim
        self.fun3dnml["boundary_output_variables"]["number_of_boundaries"] = -1
        self.fun3dnml["boundary_output_variables"]["boundary_list"] = "1-2"
        self.fun3dnml["boundary_output_variables"]["temperature"] = True
        self.fun3dnml["boundary_output_variables"]["mach"] = True
        self.fun3dnml["boundary_output_variables"]["p"] = True

        ##############################
        # fun3d settings for moving_body.input file
        self.moving_body_input = f90nml.Namelist()

        #moving body settings for funtofem to fun3d
        bodyName = self.config["csmFile"].split(".")[0]
        nBodies = 1
        nBoundaries = 1
        bndryArray = [[2]]
        bndryArray = list(bndryArray)

        #body definitions
        self.moving_body_input["body_definitions"] = f90nml.Namelist()
        self.moving_body_input["body_definitions"]["n_moving_bodies"] = nBodies
        self.moving_body_input["body_definitions"]["body_name"] = [bodyName]
        self.moving_body_input["body_definitions"]["parent_name"] = [""] # '' means motion relative to inertial ref frame
        self.moving_body_input["body_definitions"]["n_defining_bndry"] = [nBoundaries] #number of boundaries that define this body
        self.moving_body_input["body_definitions"]["defining_bndry(1,1)"] = 2 #index 1: boundary number index 2: body number
        self.moving_body_input["body_definitions"]["motion_driver"] = ["funtofem"] #tells fun3d to use motion inputs from python
        self.moving_body_input["body_definitions"]["mesh_movement"] = ["deform"] #can use 'rigid', 'deform', 'rigid+deform' with funtofem interface
        """
        pass