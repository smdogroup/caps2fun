
__all__ = ["TacsMeshLoader"]

from typing import TYPE_CHECKING
from pyNastran.bdf.bdf import read_bdf
from .material import TacsMaterial

class TacsMeshLoader:
    """
    python module to read in .bdf and/or .dat files using pyNastran
        only supports shell elements, and static loads but can create dynamic problems from it
        pynastran resources: https://pynastran-git.readthedocs.io/en/1.2/reference/bdf/cards/pyNastran.bdf.cards.html
    """
    def __init__(self, comm, dat_file:str, auto_read:bool=True, debug:bool=False):
        self._comm = comm
        self._dat_file = dat_file
        self._debug = debug
        self._bdf_info = None
        if auto_read: self.read_bdf()
    
    def read_bdf(self):
        if self._comm.rank == 0:
            debugPrint = self._debug
        else:
            debugPrint = False
        self._bdf_info = read_bdf(self._dat_file, validate=False, xref=False, debug=debugPrint)

        # make sure no subcases yet
        print(len(self.bdf_info.subcases))
        assert (len(self.bdf_info.subcases) - 1 <= 1 )

    @property
    def dat_file(self) -> str:
        return self._dat_file

    @property
    def bdf_info(self):
        assert(self._bdf_info is not None)
        return self._bdf_info

    @property
    def coordinates(self):
        return self.bdf_info.get_xyz_in_coord()

    @property
    def properties(self):
        return self.bdf_info.properties
    
    @property
    def nodes(self):
        return self.bdf_info.nodes

    @property
    def elements(self):
        return self.bdf_info.elements

    @property
    def loads(self):
        return self.bdf_info.loads

    @property
    def dvprels(self):
        return self.bdf_info.dvprels

    @property
    def materials(self):
        return self.bdf_info.materials

    @property
    def constraints(self):
        return self.bdf_info.spcs

    def process_bdf(self):

        self.process_properties()
        self.process_elements()
        self.process_nodes()
        self.process_loads()

    def process_properties(self):
        self._property_ids = []
        for idx, property_key in enumerate(self.properties):
            property = self.properties[property_key]
            #print(property)
            self._property_ids.append(property)

    def process_elements(self):
        for idx, element_key in enumerate(self.elements):
            element = self.elements[element_key]
            element_type = element.type
            property_id = element.pid
            nodes = element.nodes.copy()
            print(element_type, property_id, nodes)
        
    def process_nodes(self):
        for idx, node_key in enumerate(self.nodes):
            node = self.nodes[node_key]
            print(node.nid, node.xyz)

    def process_loads(self):
        """
        Method to process static load cards for shells
        https://pynastran-git.readthedocs.io/en/1.2/reference/bdf/cards/loads/pyNastran.bdf.cards.loads.static_loads.html
        """
        #print(self.loads)
        for load_set_key in self.loads:
            #print(f"load set idx={set_idx}")
            load_set = self.loads[load_set_key]
            #print(load_set)
            for load in load_set:
                #print(load)
                load_type = load.type
                sid = load.sid
                applied_node = load.node
                print(sid, load_type, applied_node)
                if load_type == "FORCE":
                    load_magnitude = load.mag
                    load_direction = load.xyz
                    print(load_magnitude, load_direction)
                elif load_type == "PLOAD":
                    pressure_magnitude = load.pressure 
                    print(pressure_magnitude)