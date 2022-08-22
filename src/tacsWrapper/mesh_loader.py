
__all__ = ["TacsMeshLoader"]

from typing import TYPE_CHECKING
from pyNastran.bdf.bdf import read_bdf
from mpi4py import MPI

class TacsMeshLoader:
    """
    python module to read in .bdf and .dat files using pyNastran
    """
    def __init__(self, comm, dat_file:str, debug:bool=False):
        self._comm = comm
        self._dat_file = dat_file
        self._debug = debug
        self._bdf_info = None
    
    def read_bdf(self):
        if self._comm.rank == 0:
            debugPrint = self._debug
        else:
            debugPrint = False
        self._bdf_info = read_bdf(self._dat_file, validate=False, xref=False, debug=debugPrint)
        #print(self._bdf_info) 
        self.process_bdf()


    def process_bdf(self):
        self._bdf_xpts = self._bdf_info.get_xyz_in_coord()
        print(self._bdf_xpts)

        self._bdf_properties = self._bdf_info.property_ids
        for pID in self._bdf_properties:
            print(pID)

        self.global_nastran_tacs_map()

        for tacsElementID, nastranElementID in enumerate(self._bdf_info.element_ids):
            element = self._bdf_info.elements[nastranElementID]
            elementType = element.type.upper()
            propertyID = element.pid
            #componentID = self.idMap(propertyID, self.nastranToTACSCompIDDict)
            #print(propertyID)

    def global_nastran_tacs_map(self):
        # Create Node ID map
        nastran_node_ids = self._bdf_info.node_ids
        nnodes = self._bdf_info.nnodes
        tacs_node_ids = range(nnodes)
        node_tuples = zip(nastran_node_ids, tacs_node_ids)
        self._nastran_to_tacs_node_dict = dict(node_tuples)
        print(self._nastran_to_tacs_node_dict)

        # Create Property/Component ID map
        nastran_property_ids = self._bdf_info.property_ids
        nproperties = self._bdf_info.nproperties
        tacs_property_ids = range(nproperties)
        property_tuples = zip(nastran_property_ids, tacs_property_ids)
        self._nastran_to_tacs_property_dict = dict(property_tuples)

        # Create Element ID map
        nastran_element_ids = self._bdf_info.element_ids
        nelements = self._bdf_info.nelements
        tacs_element_ids = range(nelements)
        element_tuples = zip(nastran_element_ids, tacs_element_ids)
        self._nastran_to_tacs_element_dict = dict(element_tuples)

    def idMap(self, fromIDs, tacsIDDict):
        """
        Translate fromIDs numbering from nastran numbering to tacs numbering.
        If node ID doesn't exist in nastranIDList, return -1 for entry
        """
        # Input is a list return a list
        if hasattr(fromIDs, '__iter__'):
            toIDs = [None] * len(fromIDs)
            # Iterate through list and call function recursively one element at a time
            for i, id in enumerate(fromIDs):
                toIDs[i] = self.idMap(id, tacsIDDict)
            return toIDs
        # Input is a int, return an int
        else:
            if fromIDs in tacsIDDict:
                return tacsIDDict[fromIDs]
            else:
                return -1

    def create_tacs_assembler(self):
        self.creator = tacs.TACS.Creator(self.comm, varsPerNode)