

class Function:
    """
    FUNtoFEM function base class
    """
    def __init__(self, name:str, index:int=0):
        self._name = name
        self._index = index

    @property
    def name(self) -> str:
        return self._name

    @property
    def index(self) -> int:
        return self._index

    @index.setter
    def index(self, new_index:int):
        self._index = new_index

class AerodynamicFunction:
    """
    FUNtoFEM Aerodynamic function class
    """
    def __init__(self, name:str, index:int=0):
        super(AerodynamicFunction,self).__init__(name=name, index=index)

class StructuralFunction:
    """
    FUNtoFEM Structural function class
    """
    def __init__(self, name:str, index:int=0):
        super(StructuralFunction,self).__init__(name=name, index=index)