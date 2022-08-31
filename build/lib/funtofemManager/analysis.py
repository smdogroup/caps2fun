
__all__ = ["Analysis", "Aerothermal", "Aeroelastic", "Aerothermoelastic"]

from typing import TYPE_CHECKING

class Analysis:
    """
    Base class for FUNtoFEM analysis
    """
    def __init__(self, name:str):
        self._name = name

class Aerothermal(Analysis):
    def __init__(self):
        super(Aerothermal,self).__init__(name="Aerothermal")

class Aeroelastic(Analysis):
    def __init__(self):
        super(Aeroelastic,self).__init__(name="Aeroelastic")

class Aerothermoelastic(Analysis):
    def __init__(self):
        super(Aerothermoelastic,self).__init__(name="Aerothermoelastic")