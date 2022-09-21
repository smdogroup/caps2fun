
__all__ = ["Constraint", "ZeroConstraint", "TemperatureConstraint"]


class Constraint:
    def __init__(self, 
    name:str, 
    caps_constraint:str, 
    constraint_type:str, 
    dof_constraint:int, 
    grid_displacement:float=0.0
    ):
        #assert(constraint_type in ["Displacement", "ZeroDisplacement", "Temperature"])
        dof_str = str(dof_constraint)
        for char in dof_str:
            assert(int(char) in range(0,7)) # means only allow dof 0-6, see pyMeshLoader _isDOFinString for more info
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
            "constraintType" : "Displacement",
            "dofConstraint" : self._dof_constraint,
            "gridDisplacement" : self._grid_displacement
        }

class ZeroConstraint(Constraint):
    def __init__(self, 
    name:str, 
    caps_constraint:str, 
    dof_constraint:int=123
    ):
        super(ZeroConstraint,self).__init__(name=name, caps_constraint=caps_constraint, constraint_type="Displacement", dof_constraint=dof_constraint)

    @property
    def dictionary(self) -> dict:
        return {
            "groupName" : self._caps_constraint,
            "constraintType" : "ZeroDisplacement",
            "dofConstraint" : self._dof_constraint
        }

class TemperatureConstraint(Constraint):
    def __init__(self,
    name:str,
    caps_constraint:str,
    temperature:float=300.0
    ):
        super(TemperatureConstraint,self).__init__(
            name=name,
            caps_constraint=caps_constraint,
            constraint_type="Thermal",
            dof_constraint=0,
            grid_displacement=temperature
        )

    @property
    def dictionary(self) -> dict:
        return super(TemperatureConstraint,self).dictionary