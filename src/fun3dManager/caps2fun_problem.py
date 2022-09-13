

__all__ = ["CapsFun3d"]

from typing import TYPE_CHECKING
from capsManager.pointwise.pointwise_aim import PointwiseAim
from capsManager.fun3d.fun3d_aim import Fun3dAim
import os

class CapsFun3d:
    def __init__(self, pointwise_aim:PointwiseAim, fun3d_aim:Fun3dAim):
        assert(pointwise_aim.is_setup)
        assert(fun3d_aim.is_setup)
        self._pointwise_aim = pointwise_aim 
        self._fun3d_aim = fun3d_aim

        # add fun3d BC identical to pointwise
        self.fun3d_aim.set_boundary_condition(self.pointwise_aim.boundary_condition)

    @property
    def pointwise_aim(self) -> PointwiseAim:
        return self._pointwise_aim

    @property
    def fun3d_aim(self) -> Fun3dAim:
        return self._fun3d_aim

    def build_mesh(self):
        """
        call pointwise mesher to build the mesh.ugrid file
        """
        self.pointwise_aim.run_pointwise()
        # could have alternate mesh builders like tetgen following this

        # link the fluid mesh from pointwise to fun3d
        self.fun3d_aim.aim.input["Mesh"].link(self.pointwise_aim.aim.output["Volume_Mesh"])

    def prepare_fun3d(self, add_to_map_bc:bool=False):
        """
        prepare fun3d by building 
        """

        # run the fun3d pre analysis to build namelist and mapbc files
        self.fun3d_aim.pre_analysis()

        # add appropriate bc names to mapbc file
        if add_to_map_bc:
            self._write_map_bc()

        self.fun3d_aim.write()

        """
        for complex step have to add perturb.input file
        """
        # #get the caps2fun project dir
        # caps2fun_proj_dir = os.environ["CAPS2FUN"]
        # archive_folder = os.path.join(caps2fun_proj_dir, "archive")

        # #move perturb.input from archive folder if complex mode
        # if (self.complex):
        #     src = os.path.join(archive_folder,"perturb.input")
        #     dest = os.path.join(caps_flow_dir, "perturb.input")
        #     shutil.copy(src, dest)

    def _write_map_bc(self):
        #add names of the BCs to the mapbc file
        mapbc = os.path.join(self.fun3d_aim.flow_directory, "fun3d_CAPS.mapbc")
        mapbc_hdl = open(mapbc, "r")
        lines = mapbc_hdl.readlines()
        mapbc_hdl.close()
        wr_hdl = open(mapbc,"w")
        # TODO : make these settings visible from csm file
        bc_inds = ["1","2"]
        bc_names = ["Farfield", "wall"]
        ind = 0
        for line in lines:
            chunks = line.split(" ")
            if (len(chunks) > 1):
                line = line.strip()
                if (bc_inds[ind] in chunks[0]):
                    line += " " + bc_names[ind]
                    ind += 1
                line += "\n"

            wr_hdl.write(line)
        wr_hdl.close()


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
