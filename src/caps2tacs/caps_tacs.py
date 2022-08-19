
__all__ = ["CapsTacs"]

from capsWrapper.tacs import TacsAim
from capsWrapper.egads import EgadsAim
import os


class CapsTacs:
    """
    Module to handle caps and tacs interface problems
    """
    def __init__(self, tacs_aim : TacsAim, egads_aim : EgadsAim):
        """
        Module to handle caps and tacs interface problems
            tacs_aim : provide a fully setup tacs aim wrapper object
            egads_aim : provide a fully setup egads aim wrapper object
        """

        # make sure the aim wrapper modules are setup properly before proceeding
        assert(tacs_aim.is_setup)
        assert(egads_aim.is_setup)

        # get the aims from each of my wrapper aim classes
        self._tacs_aim = tacs_aim.aim
        self._egads_aim = egads_aim.aim

        self._dat_file = None

        # link meshes and run pre analysis
        self.link_meshes()
        self.tacs_pre_analysis()

    @property
    def tacs_aim(self):
        return self._tacs_aim

    @property
    def egads_aim(self):
        return self._egads_aim

    @property
    def dat_file(self) -> str:
        return self._dat_file

    def link_meshes(self):
        # link the meshes
        self.tacs_aim.input["Mesh"].link(self.egads_aim.output["Surface_Mesh"])

    def tacs_pre_analysis(self):
        # run the tacs aim preanalysis
        self.tacs_aim.preAnalysis()

        self._dat_file = os.path.join(self.tacs_aim.analysisDir, self.tacs_aim.input.Proj_Name + '.dat')

