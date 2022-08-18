

__all__ = ["TacsAim"]

from typing import TYPE_CHECKING
import pyCAPS
from .material import Material

class TacsAim:
    """
    Wrapper class for TacsAim with default build setting in different scenarios
    applies default settings and spits it back out at the end
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

    @property
    def is_built(self) -> bool:
        return self._file_format and self._analysis and self._materials and self._constraints

    def add_material(self, material:Material):
        """
        add a Material object to the materials list
        """
        self._materials.append(material)

    def setup_aim(self, large_format:bool=True, static:bool=True, ):
        
        #increase the precision in the BDF file
        self._aim.input.File_Format = "Large" if large_format else "Small"
        self._aim.input.Mesh_File_Format = "Large" if large_format else "Small"

        self._aim.input.Analysis_Type = "Static" if static else "Dynamic"

        # add materials to tacsAim
        self._aim.input.Material = {material.name : material.dictionary for material in self._materials}

        # note that setup is finished now
        self._setup = True
        
    def set_property(self):
        # Material properties section
        shell  = {"propertyType" : "Shell",
                "membraneThickness" : 0.006,
                "material"        : "madeupium",
                "bendingInertiaRatio" : 1.0, # Default
                "shearMembraneRatio"  : 5.0/6.0} # Default
        
        self._aim.input.Property = {"plate": shell}
        
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
        returns the pre-built tacsAim object
        """
        return self._aim