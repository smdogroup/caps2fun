
__all__ = ["PointwiseAim"]

from typing import TYPE_CHECKING
import pyCAPS

class PointwiseAim:
    def __init__(self, caps_problem:pyCAPS.Problem):
        self._aim = None