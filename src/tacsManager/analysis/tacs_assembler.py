
__all__ = ["TacsAssembler"]

from tacsManager.nastran_reader.mesh_loader import TacsMeshLoader
from tacs.TACS import Creator

class TacsAssembler:
    """
    module to make an initial Tacs Assembler
    """
    def __init__(self, mesh_loader:TacsMeshLoader, comm, vars_per_node:int=6):
        self._mesh_loader = mesh_loader
        self.comm = comm
        self.vars_per_node = vars_per_node

        self._creator = None
    
    @property
    def mesh_loader(self) -> TacsMeshLoader:
        return self._mesh_loader

    @property
    def assembler(self):
        """
        method to get tacs assembler
        """
        self._creator = Creator(self.comm, self.vars_per_node)
        if self.comm.rank == 0:
            print(self._creator)


