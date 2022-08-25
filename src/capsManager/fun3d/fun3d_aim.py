
__all__ = ["Fun3dAim"]


from typing import TYPE_CHECKING
import pyCAPS

class Fun3dAim:
    def __init__(self, caps_problem:pyCAPS.Problem):
        self._aim = None