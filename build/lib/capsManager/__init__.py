"""
Manager package for ESP/CAPS AIMs (Analyses Interface Modules) such as tacsAim, egadsTessAim, pointwiseAim, fun3dAim
a caps problem is initialized with CapsProblem module using a csm file provided by the user
"""
from .egads import *
from capsManager.file_manager import *
from capsManager.fun3d import *
from capsManager.pointwise import *
from capsManager.caps_problem import *
from capsManager.tacs import *
