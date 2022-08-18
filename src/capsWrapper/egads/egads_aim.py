

__all__ = ["EgadsAim"]

import pyCAPS

class EgadsAim:
    def __init__(self, caps_problem:pyCAPS.Problem):
        self._aim = caps_problem.analysis.create(aim="egadsTessAIM")