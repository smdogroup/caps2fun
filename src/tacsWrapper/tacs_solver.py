
__all__ = ["TacsSolver"]

from typing import TYPE_CHECKING
from pyNastran.bdf.bdf import read_bdf
from mpi4py import MPI

# TODO : write alternate tacs wrapper to pytacs
class TacsSolver:
    """
    python interface to TACS, aka alternate TACS Wrapper class to pytacs
    """
    def __init__(self, dat_file:str):
        self._dat_file = dat_file