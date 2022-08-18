
__all__ = ["Parameter", "ConfigParameter", "DesignParameter"]

from csmWriter.statements.base_statement import Statement

class Parameter:
    def __init__(self, name:str, value:float=0.0):
        self._name = name
        self._value = value
    
    @property
    def name(self) -> str:
        return self._name

    @property
    def value(self) -> float:
        return self._value

    def __str__(self):
        return self.name

class ConfigParameter(Parameter):
    def __init__(self, name:str, value:int):
        super(ConfigParameter, self).__init__(name=name, value=value)
    
    @property
    def statement(self) -> Statement:
        return Statement(content=f"cfgpmtr\t{self.name}\t{self.value}")

class DesignParameter(Parameter):
    def __init__(self, name:str, value:float):
        super(DesignParameter, self).__init__(name=name, value=value)

    @property
    def statement(self) -> Statement:
        return Statement(content=f"despmtr\t{self.name}\t{self.value}")

    