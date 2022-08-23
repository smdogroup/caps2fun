
__all__ = ["Load", "Pressure", "GridForce"]

from typing import TYPE_CHECKING, List

class Load:
    """
    generates load dictionary for tacs Aim
    load = {"groupName" : "plate",
            "loadType" : "Pressure",
            "pressureForce" : pressload}
    caps_load
        matches caps load attribute in CSM file
    """
    def __init__(self, name:str, caps_load:str, load_type:str):
        assert(load_type in ["GridForce", "GridMoment", "Rotational", "Thermal", "Pressure", "PressureDistribute", "PressureExternal", "Gravity"])
        self._name = name
        self._caps_load = caps_load
        self._load_type = load_type

    @property
    def name(self) -> str:
        return self._name

class Pressure(Load):
    def __init__(self, name:str, caps_load:str, force:float=2.0E6):
        super(Pressure,self).__init__(name=name, caps_load=caps_load, load_type="Pressure")
        self._pressure_force = force

    @property
    def dictionary(self) -> dict:
        return {
            "groupName" : self._caps_load,
            "loadType" : self._load_type,
            "pressureForce" : self._pressure_force
        }

class GridForce(Load):
    def __init__(self, name:str, caps_load:str, direction:List[float]=[0.,0.,1.], magnitude:float=1.0E3):
        super(GridForce,self).__init__(name=name, caps_load=caps_load, load_type="GridForce")
        self._direction_vector = direction
        self._force_scale_factor = magnitude

    @property
    def dictionary(self) -> dict:
        return {
            "groupName" : self._caps_load,
            "loadType" : self._load_type,
            "directionVector" : self._direction_vector,
            "forceScaleFactor" : self._force_scale_factor
        }

    