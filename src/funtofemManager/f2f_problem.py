
__all__ = ["F2FProblem"]

from typing import TYPE_CHECKING

class F2FProblem:
    def __init__(self, name:str, id:int=0):
        """
        FUNtoFEM problem class, holds all the data for coupled FUNtoFEM simulation
            To create a model, instantiate it, and then add bodies and scenarios to it.
            The functions and derivative values can be extracted from the FUNtoFEM problem object
        """
        self._name = name
        self._id = id

        self._scenarios = []
        self._bodies = []

    def add_body(self, new_body):
        self._bodies.append(new_body)