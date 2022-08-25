
import f90nml
__all__ = ["Fun3dAim"]


from typing import TYPE_CHECKING
import pyCAPS

class Fun3dAim:
    def __init__(self, caps_problem:pyCAPS.Problem):
        self._aim = None

    # self.fun3dAim.input.Boundary_Condition = {"wall": self.wallBC,
    #             "Farfield": {"bcType":"Farfield"}}

    #     #add thickDVs and geomDVs to caps
    #     DVdict = {}
    #     for DV in self.DVdict:
    #         dvname = DV["name"]
    #         if (DV["type"] == "shape"): #geomDV, add empty entry into DV dicts
    #             DVdict[dvname] = {}

    #     #input DVdict and DVRdict into tacsAim
    #     if (len(DVdict) > 0): self.fun3dAim.input.Design_Variable = DVdict

    #     #fun3d design sensitivities settings
    #     self.fun3dAim.input.Design_SensFile = True
    #     self.fun3dAim.input.Design_Sensitivity = True

    #     #############################
        