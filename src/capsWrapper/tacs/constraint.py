
__all__ = ["Constraint", "ZeroConstraint"]

from doctest import DONT_ACCEPT_TRUE_FOR_1


class Constraint:
    def __init__(self, name:str, caps_constraint:str, constraint_type:str, dof_constraint:int, grid_displacement:float=0.0):
        assert(constraint_type in ["Displacement", "ZeroDisplacement"])
        dof_str = str(dof_constraint)
        for char in dof_str:
            assert(int(char) in range(1,7)) # means only allow dof 1,2,3,4,5,6
        self._name = name
        self._caps_constraint = caps_constraint
        self._dof_constraint = dof_constraint
        self._grid_displacement = grid_displacement

    @property
    def name(self) -> str:
        return self._name

    @property
    def dictionary(self) -> dict:
        return {
            "groupName" : self._caps_constraint,
            "dofConstraint" : self._dof_constraint,
            "gridDisplacement" : self._grid_displacement
        }

class ZeroConstraint(Constraint):
    def __init__(self, name:str, caps_constraint:str, dof_constraint:int=123):
        super(ZeroConstraint,self).__init__(name=name, caps_constraint=caps_constraint, constraint_type="Displacement", dof_constraint=dof_constraint)

    @property
    def dictionary(self) -> dict:
        return {
            "groupName" : self._caps_constraint,
            "dofConstraint" : self._dof_constraint
        }