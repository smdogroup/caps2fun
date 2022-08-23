
__all__ = ["Property", "ShellProperty"]

from .material import Material
from typing import TYPE_CHECKING


class Property:
    def __init__(self, caps_group:str, material:Material, property_type:str):
        self._caps_group = caps_group
        self._material = material
        self._property_type = property_type

    @property
    def capsGroup(self) -> str:
        """
        return capsGroup attribute associated with this property
        """
        return self._caps_group

    @property
    def dictionary(self) -> dict:
        """
        return property dictionary, however this is only fully defined in subclasses
        """
        return {}

class ShellProperty(Property):
    """
    Example ShellProperty Dictionary
    shell  = {"propertyType" : "Shell",
                    "membraneThickness" : 0.006,
                    "material"        : "madeupium",
                    "bendingInertiaRatio" : 1.0, # Default
                    "shearMembraneRatio"  : 5.0/6.0} # Default

    self._aim.input.Property = {"plate": shell}
    """
    # TODO : add other available settings for shell properties -> mass, inertias, etc.
    def __init__(self, caps_group:str, material:Material, membrane_thickness:float, bending_inertia:float=1.0, shear_membrane_ratio:float=5.0/6.0):
        super(ShellProperty,self).__init__(caps_group=caps_group, material=material, property_type="Shell")
        self._membrane_thickness = membrane_thickness
        self._bending_inertia = bending_inertia
        self._shear_membrane_ratio = shear_membrane_ratio

    @property
    def membrane_thickness(self) -> float:
        return self._membrane_thickness

    @membrane_thickness.setter
    def membrane_thickness(self, new_thickness:float):
        self._membrane_thickness = new_thickness

    @property
    def dictionary(self) -> dict:
        """
        return property dictionary to pass into tacsAim
        """
        return {"propertyType" : self._property_type,
                "membraneThickness" : self._membrane_thickness,
                "material"        : self._material.name,
                "bendingInertiaRatio" : self._bending_inertia, # Default
                "shearMembraneRatio"  : self._shear_membrane_ratio} # Default


