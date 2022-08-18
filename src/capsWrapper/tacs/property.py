
from .material import Material

"""
shell  = {"propertyType" : "Shell",
                "membraneThickness" : 0.006,
                "material"        : "madeupium",
                "bendingInertiaRatio" : 1.0, # Default
                "shearMembraneRatio"  : 5.0/6.0} # Default

self._aim.input.Property = {"plate": shell}
"""

class Property:
    def __init__(self, name:str, material:Material, property_type:str):
        self._name = name
        self._material = material
        self._property_type = property_type

    @property
    def name(self) -> str:
        return self._name

class ShellProperty:
