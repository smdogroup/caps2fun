
__all__ = ["Variable", "StructuralVariable", "AerodynamicVariable"]

from typing import TYPE_CHECKING

class Variable:
    def __init__(self, name:str, value:float=0.0, lower_bound:float=0.0, upper_bound:float=0.0, scale:float=1.0,
        active:bool=True, coupled:bool=False, id:int=0):
        """
        FUNtoFEM base variable class
        """
        self._name = name
        self._value = value
        self._lower_bound = lower_bound
        self._upper_bound = upper_bound
        self._scale = scale
        self._active = active
        self._coupled = coupled
        self._id = id

    @property
    def name(self) -> str:
        return self._name

    @property
    def value(self) -> float:
        return self._value

    @value.setter
    def value(self, new_value:float):
        self._value = new_value

    @property
    def active(self) -> bool:
        return self._active

    @property
    def coupled(self) -> bool:
        return self._coupled

    @property
    def id(self) -> int:
        return self._id

class AerodynamicVariable(Variable):
    def __init__(self, name:str, value:float=0.0, lower_bound:float=0.0, upper_bound:float=0.0, scale:float=1.0,
        active:bool=True, coupled:bool=False, id:int=0):
        """
        FUNtoFEM Aerodynamic variable class
        """
        super(AerodynamicVariable,self).__init__(
            name=name, 
            value=value, 
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            scale=scale,
            active=active,
            coupled=coupled,
            id=id
            )

class StructuralVariable(Variable):
    def __init__(self, name:str, value:float=0.0, lower_bound:float=0.0, upper_bound:float=0.0, scale:float=1.0,
        active:bool=True, coupled:bool=False, id:int=0):
        """
        FUNtoFEM Structural variable class
        """
        super(StructuralVariable,self).__init__(
            name=name, 
            value=value, 
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            scale=scale,
            active=active,
            coupled=coupled,
            id=id
            )