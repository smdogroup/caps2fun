"""
Manager package for ESP/CAPS AIMs (Analyses Interface Modules) such as tacsAim, egadsTessAim, pointwiseAim, fun3dAim
a caps problem is initialized with CapsProblem module using a csm file provided by the user
"""
from .egads import *
from .file_manager import *
from .fun3d import *
from .pointwise import *
from .caps_problem import *
from .tacs import *