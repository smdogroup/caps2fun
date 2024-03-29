"""
Sean Engelstad
Georgia Tech SMDO
08/18/2022
"""

__all__ = ["CapsProblem", "CapsStruct", "CapsFluid"]

import pyCAPS
from typing import TYPE_CHECKING

# import each of the aim modules here
from capsManager.tacs.tacs_aim import TacsAim
from capsManager.egads.egads_aim import EgadsAim
from capsManager.pointwise.pointwise_aim import PointwiseAim
from capsManager.fun3d.fun3d_aim import Fun3dAim

class CapsProblem:
    """
    Wrapper class to make a pycaps problem and export aims
    """
    def __init__(self, problem:pyCAPS.Problem):
        self._caps_problem = problem

    @classmethod
    def default(cls, csmFile:str, problemName:str="capsDir"):
        """
        auto build a caps problem
        """
        problem = pyCAPS.Problem(problemName=problemName, capsFile=csmFile, outLevel=1)
        return cls(problem)

    @property
    def geometry(self):
        return self._caps_problem.geometry

    @property
    def view(self):
        self.geometry.view()

    @property
    def design_parameters(self):
        return self.geometry.despmtr.keys()

class CapsStruct(CapsProblem):
    """
    Base class for Structure problems with ESP/CAPS
    Often uses TACS for structure solver
    """
    def __init__(self, problem:pyCAPS.Problem):
        super(CapsStruct,self).__init__(problem=problem)

    @classmethod
    def default(cls, csmFile:str, problemName:str="capsStruct"):
        """
        auto build a caps struct problem
        """
        problem = pyCAPS.Problem(problemName=problemName, capsFile=csmFile, outLevel=1)
        return cls(problem)

    @property
    def tacsAim(self) -> TacsAim:
        return TacsAim(caps_problem=self._caps_problem)

    @property
    def egadsAim(self) -> EgadsAim:
        return EgadsAim(caps_problem=self._caps_problem)

    @property
    def design_parameters(self):
        return self.geometry.despmtr.keys()

class CapsFluid(CapsProblem):
    """
    Base class for Aerodynamics/Fluid problems with ESP/CAPS
    Often uses Fun3d for solver and pointwise for meshing
    """
    def __init__(self, problem:pyCAPS.Problem):
        super(CapsFluid,self).__init__(problem=problem)

        # set to not overwrite fun3d nml analysis
        self._caps_problem.analysis["fun3d"].input.Overwrite_NML = False

    @classmethod
    def default(cls, csmFile:str, problemName:str="capsFluid"):
        """
        auto build a caps struct problem
        """
        problem = pyCAPS.Problem(problemName=problemName, capsFile=csmFile, outLevel=1)
        return cls(problem)

    @property
    def pointwiseAim(self) -> PointwiseAim:
        return PointwiseAim(caps_problem=self._caps_problem)

    @property
    def fun3dAim(self) -> Fun3dAim:
        return Fun3dAim(caps_problem=self._caps_problem)
