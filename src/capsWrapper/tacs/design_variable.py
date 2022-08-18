
__all__ = ["ShapeVariable", "ThicknessVariable"]

from .property import ShellProperty
from .material import Material

class ShapeVariable:
    """
    shape variables are ESP/CAPS despmtr variables, etc.
    """
    def __init__(self, name:str):
        self._name = name

    @property
    def name(self) -> str:
        return self._name
    
    @property
    def DV_dictionary(self) -> dict:
        return {}

    @property
    def is_thickness_DV(self) -> bool:
        return False

class ThicknessVariable:
    """
    shape variables are ESP/CAPS despmtr variables, etc.
    """
    def __init__(self, name:str, caps_group:str, value:float, material:Material=None):
        self._name = name
        self._caps_group = caps_group
        self._value = value
        self._material = material

    @property
    def name(self) -> str:
        return self._name
    
    @property
    def DV_dictionary(self) -> dict:
        """
        Design Variable Dictionary
        """
        return {"groupName" : self._caps_group,
                "initialValue" : self._value,
                "lowerBound" : self._value*0.5,
                "upperBound" : self._value*1.5,
                "maxDelta"   : self._value*0.1} 

    @property
    def DVR_dictionary(self) -> dict:
        """
        Design Variable Relation Dictionary
        """
        return {"variableType": "Property",
                "fieldName" : "T",
                "constantCoeff" : 0.0,
                "groupName" : self._name,
                "linearCoeff" : 1.0}

    @property
    def has_material(self) -> bool:
        return self._material is not None

    @property
    def shell_property(self) -> ShellProperty:
        assert(self._material is not None)
        return ShellProperty(caps_group=self._caps_group, material=self._material, membrane_thickness=self._value)

    @property
    def is_thickness_DV(self) -> bool:
        return True
