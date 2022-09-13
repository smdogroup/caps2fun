

import numpy as np
import os, zmq
from mpi4py import MPI
from funtofem import TransferScheme
from funtofemManager.f2f_problem import Problem
from fun3d.mda.fsa.fun3d_aero import (Fun3dAero, Fun3dAeroException,
                                      SurfaceMesh, AeroLoads, AdjointProduct)

from typing import TYPE_CHECKING

class Fun3dClient:
    """
    FUNtoFEM client class for FUN3D. Works for both steady and unsteady analysis.
    Requires the FUN3D directory structure.
        For the forward analysis, the FUN3D interface will operate in the scenario.name/Flow directory.
        For the adjoint analysis, the FUN3D interface will operate in the scenario.name/Adjoint directory.
    
    To tell FUN3D that a body's motion should be driven by FUNtoFEM, set *motion_driver(i)='funtofem'
    """
    def __init__(self, comm, problem:Problem, flow_dt:float=1.0,
                 host:float="localhost", port_base:int=49200):

        self._comm = comm
        self._flow_dt = flow_dt

        # instantiate fun3d
        mpi_rank = MPI.COMM_WORLD.Get_rank()
        port = port_base + mpi_rank
        context = zmq.Context()
        self._fun3d_client = Fun3dAero.Client(
            context,
            f"tcp://{host}:{port}",
            zmq.REQ
        )

        # get initial aerodynamic and surface meshes