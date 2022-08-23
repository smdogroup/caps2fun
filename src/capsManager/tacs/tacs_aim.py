

__all__ = ["TacsAim"]

from typing import TYPE_CHECKING, List
import pyCAPS
from .material import Material
from .property import ShellProperty
from .constraint import Constraint
from .load import *
from .design_variable import ShapeVariable, ThicknessVariable

class TacsAim:
    """
    Wrapper class for TacsAim with default build setting in different scenarios
    applies default settings and spits it back out at the end
    only supports shell properties at the moment
    """
    
    def __init__(self, caps_problem:pyCAPS.Problem):
        self._aim = caps_problem.analysis.create(aim = "tacsAIM", name = "tacs")

        self._materials = []
        self._loads = []
        self._properties = []
        self._constraints = []
        self._design_variables = []

        # build flags
        self._setup = False

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

    def setup_aim(self, large_format:bool=True, static:bool=True):
        
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

        # add the design variables to the DesignVariable and DesignVariableRelation properties
        self._aim.input.Design_Variable_Relation = {dv.name : dv.DVR_dictionary for dv in self._design_variables if isinstance(dv, ThicknessVariable)}
        self._aim.input.Design_Variable = {dv.name : dv.DV_dictionary for dv in self._design_variables}

        # note that setup is finished now
        self._setup = True

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