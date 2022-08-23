

__all__ = ["TacsAim"]

from typing import TYPE_CHECKING, List
import pyCAPS
import os
from capsManager.tacs.material import Material
from capsManager.tacs.property import ShellProperty
from capsManager.tacs.constraint import Constraint
from capsManager.tacs.load import *
from capsManager.tacs.design_variable import ShapeVariable, ThicknessVariable

class TacsAim:
    """
    Wrapper class for TacsAim with default build setting in different scenarios
    applies default settings and spits it back out at the end
    only supports shell properties at the moment
    """
    
    def __init__(self, caps_problem:pyCAPS.Problem):
        self._aim = caps_problem.analysis.create(aim = "tacsAIM", name = "tacs")
        
        # geometry and design parameters to change the design of the CSM file during an optimization
        self._design_parameters = caps_problem.geometry.despmtr.keys()
        self._geometry = caps_problem.geometry

        self._materials = []
        self._loads = []
        self._properties = []
        self._constraints = []
        self._design_variables = []

        # build flags
        self._setup = False
        self._first_setup = True

    def update_design(self, design_dict:dict):
        """
        method to change the values of each design variable in tacs, caps
        input x is a dictionary of values for each variable {"name" : value}
        """
        # change setup to False since we have a new design we have to set up the aim again
        self._setup = False

        for design_variable in self._design_variables:
            if design_variable.name in design_dict:
                new_value = design_dict[design_variable.name]
                if isinstance(design_variable, ShapeVariable):
                    self._geometry.despmtr[new_value].value = new_value
                elif isinstance(design_variable, ThicknessVariable):
                    design_variable.value = new_value
                    # change that property
                    matching_property = False
                    for property in self._properties:
                        if isinstance(property, ShellProperty):
                            matching_property = property.membrane_thickness == design_variable.caps_group
                            if matching_property:
                                property.membrane_thickness = new_value
                    if not matching_property:
                        raise AssertionError("Couldn't find matching property...")

    def add_material(self, material:Material):
        """
        add a Material object to the materials list
        """
        assert(isinstance(material, Material))
        self._materials.append(material)

    def add_variable(self, variable):
        """
        automatically adds Shape and/or Thickness Variables
        """
        self._design_variables.append(variable)
        if isinstance(variable, ThicknessVariable):
            if variable.has_material:
                self._properties.append(variable.shell_property) 
            else:
                raise AssertionError("Did not set material for this property...")       

    def add_property(self, property:ShellProperty):
        assert(isinstance(property, ShellProperty))
        self._properties.append(property)        

    def add_constraint(self, constraint:Constraint):
        assert(isinstance(constraint, Constraint))
        self._constraints.append(constraint)    

    def add_load(self, load:Load):
        assert(isinstance(load, Load))
        self._loads.append(load)   

    def setup_aim(self, large_format:bool=True, static:bool=True, auto_shape_variables:bool=True):
        
        # make sure there is at least one material, property, constraint, etc.
        assert(len(self._materials) > 0)
        assert(len(self._properties) > 0)
        assert(len(self._constraints) > 0)

        #increase the precision in the BDF file
        self._aim.input.File_Format = "Large" if large_format else "Small"
        self._aim.input.Mesh_File_Format = "Large" if large_format else "Small"

        # set the analysis type
        if static:
            self._aim.input.Analysis_Type = "Static"
        else:
            raise AssertionError("Analysis types other than static analyses for tacsAim are not supported yet.")

        # add materials to tacsAim
        self._aim.input.Material = {material.name : material.dictionary for material in self._materials}

        # add properties to tacsAim
        self._aim.input.Property = {prop.capsGroup : prop.dictionary for prop in self._properties}

        # add constraints to tacsAim
        self._aim.input.Constraint = {con.name : con.dictionary for con in self._constraints}

        # add loads to tacsAim
        if len(self._loads) > 0:
            self._aim.input.Load = {load.name : load.dictionary for load in self._loads}

        if auto_shape_variables and self._first_setup:
            for despmtr in self._design_parameters:
                shape_var = ShapeVariable(name=despmtr)
                self.add_variable(variable=shape_var)
            self._first_setup = False

        # add the design variables to the DesignVariable and DesignVariableRelation properties
        self._aim.input.Design_Variable_Relation = {dv.name : dv.DVR_dictionary for dv in self._design_variables if isinstance(dv, ThicknessVariable)}
        self._aim.input.Design_Variable = {dv.name : dv.DV_dictionary for dv in self._design_variables}

        # note that setup is finished now
        self._setup = True

    @property
    def design_variables(self) -> List[ShapeVariable or ThicknessVariable]:
        return self._design_variables

    def pre_analysis(self):
        """
        provide access to the tacs aim preAnalysis for running 
        """
        self.aim.preAnalysis()

    @property
    def analysis_dir(self) -> str:
        return self.aim.analysisDir

    @property
    def dat_file(self) -> str:
        return self.project_name + '.dat'

    @property
    def dat_file_path(self) -> str:
        return os.path.join(self.analysis_dir, self.dat_file)

    @property
    def sens_file(self) -> str:
        return self.project_name + ".sens"

    @property
    def sens_file_path(self) -> str:
        return os.path.join(self.analysis_dir, self.sens_file)

    @property
    def project_name(self) -> str:
        return self.aim.input.Proj_Name

    def post_analysis(self):
        self.aim.postAnalysis()

    @property
    def is_setup(self) -> bool:
        return self._setup 
        
    def set_constraint(self):
        # constraint section
        constraint = {"groupName" : "edge",
                    "dofConstraint" : 123456}
        
        self._aim.input.Constraint = {"edgeConstraint": constraint}
    
    def set_load(self):
        #loads section
        pressload = 1.0e5
        
        # Set load
        load = {"groupName" : "plate",
                "loadType" : "Pressure",
                "pressureForce" : pressload}
        
        # Set loads
        self._aim.input.Load = {"appliedPressure": load }

    @property
    def aim(self):
        """
        returns the auto-built tacsAim object
        """
        return self._aim