
__all__ = ["Derivative", "Gradient"]

from typing import TYPE_CHECKING, List
from funtofemManager.function import Function
from funtofemManager.variable import Variable, StructuralVariable, AerodynamicVariable

class Derivative:
    def __init__(self, function:Function, variable:Variable):
        self._function = function
        self._variable = variable

        self._value = None

    @property
    def value(self) -> float:
        return self._value

    @value.setter
    def value(self, new_value:float):
        self._value = new_value

    @property
    def function(self) -> Function:
        return self._function

    @property
    def function_name(self) -> str:
        return self._function.name

    @property
    def variable(self) -> Variable:
        return self._variable

    @property
    def variable_name(self) -> str:
        return self.variable.name

    @property
    def aerodynamic(self) -> bool:
        return isinstance(self.variable, AerodynamicVariable)

    @property
    def structural(self) -> bool:
        return isinstance(self.variable, StructuralVariable)

class Gradient:
    """
    FUNtoFEM gradient class, holds a group of derivative objects for each function
    List of gradients, one for each function will represent 
    """
    def __init__(self, function:Function, derivatives:List[Derivative]):

        # ensure the derivatives all are for the same function object
        for derivative in derivatives:
            assert(derivative.function == function)

        self._function = function
        self._derivatives = derivatives

    @property
    def function(self) -> Function:
        return self._function

    @property
    def derivatives(self) -> List[Derivative]:
        return self._derivatives

    @property
    def aerodynamic_derivatives(self) -> List[Derivative]:
        return [deriv for deriv in self._derivatives if isinstance(deriv.variable, AerodynamicVariable)]

    @property
    def structural_derivatives(self) -> List[Derivative]:
        return [deriv for deriv in self._derivatives if isinstance(deriv.variable, StructuralVariable)]

