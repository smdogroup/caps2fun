
__all__ = ["TacsSolver"]

from typing import TYPE_CHECKING
from .nastran_reader.mesh_loader import TacsMeshLoader
from .analysis.tacs_assembler import TacsAssembler
from mpi4py import MPI

# TODO : write alternate tacs wrapper to pytacs
class TacsSolver:
    """
    python interface to TACS, aka alternate TACS Wrapper class to pytacs
    """
    def __init__(self, comm, dat_file:str):
        self._dat_file = dat_file
        self._mesh_loader = TacsMeshLoader(comm=comm, dat_file=dat_file)
        self._assembler = TacsAssembler(mesh_loader=self._mesh_loader, comm=comm).assembler

    # def create_static_problem(self):
    #     assert(self._assembler is not None)
    #     problem = self.createStaticProblem(name)

    #     if self.assembler is None:
    #         raise self._initializeError()

    #     problem = tacs.problems.static.StaticProblem(name, self.assembler, self.comm,
    #                                                  self.outputViewer, self.meshLoader, options)
    #     # Set with original design vars and coordinates, in case they have changed
    #     problem.setDesignVars(self.x0)
    #     problem.setNodes(self.Xpts0)

    #     if 'LOAD' in subCase.params:
    #         loadsID = subCase.params['LOAD'][0]
    #         # Get loads and scalers for this load case ID
    #         loadSet, loadScale, _ = self.bdfInfo.get_reduced_loads(loadsID)
    #         # Loop through every load in set and add it to problem
    #         for loadInfo, scale in zip(loadSet, loadScale):
    #             # Add any point force or moment cards
    #             if loadInfo.type == 'FORCE' or loadInfo.type == 'MOMENT':
    #                 nodeID = loadInfo.node_ref.nid

    #                 loadArray = numpy.zeros(vpn)
    #                 if loadInfo.type == 'FORCE' and vpn >= 3:
    #                     loadArray[:3] += scale * loadInfo.scaled_vector
    #                 elif loadInfo.type == 'MOMENT' and vpn >= 6:
    #                     loadArray[3:6] += scale * loadInfo.scaled_vector
    #                 problem.addLoadToNodes(nodeID, loadArray, nastranOrdering=True)

    #             # Add any gravity loads
    #             elif loadInfo.type == 'GRAV':
    #                 inertiaVec = np.zeros(3, dtype=self.dtype)
    #                 inertiaVec[:3] = scale * loadInfo.scale * loadInfo.N
    #                 problem.addInertialLoad(inertiaVec)

    #             # Add any pressure loads
    #             # Pressure load card specific to shell elements
    #             elif loadInfo.type == 'PLOAD2':
    #                 elemIDs = loadInfo.eids
    #                 pressure = scale * loadInfo.pressure
    #                 problem.addPressureToElements(elemIDs, pressure, nastranOrdering=True)

    #             # Alternate more general pressure load type
    #             elif loadInfo.type == 'PLOAD4':
    #                 self._addPressureFromPLOAD4(problem, loadInfo, scale)

    #             else:
    #                 self._TACSWarning("Unsupported load type "
    #                             f" '{loadInfo.type}' specified for load set number {loadInfo.sid}, skipping load")

    #         # append to list of structural problems
    #         structProblems[subCase.id] = problem