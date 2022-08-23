
__all__ = ["PytacsFunction", "MassStress"]

from typing import TYPE_CHECKING, List
from tacs import functions, pyTACS
from abc import ABC, abstractmethod
import os

class PytacsFunction(ABC):
    """
    Base class for pytacs analysis functions
        can create your own custom analysis functions as well and use them in caps2tacs
        however need it to be subclass so that it will provide these properties -> coordinate derivatives, etc.
    """
    def __init__(self, name:str):
        self._name = name
        self._fea_solver = None
        self._function_names = []
        self._analysis_dir = None
        self._write_idx = 0

    @property
    def name(self) -> str:
        return self._name

    @property
    def f5_filename(self) -> str:
        """
        name of f5 file without .f5 extension
        """
        return f"{self.name}_{self._write_idx}"

    @property
    def f5_filepath(self) -> str:
        load_set_str = "000"
        return os.path.join(self.analysis_dir, f"{self.f5_filename}_{load_set_str}.f5")

    @property
    def num_functions(self) -> int:
        return len(self.function_names)

    @property
    def fea_solver(self):
        """
        pytacs fea solver object
        """
        assert(self._fea_solver is not None)
        return self._fea_solver

    @fea_solver.setter
    def fea_solver(self, new_solver):
        self._fea_solver = new_solver

    @property
    def num_nodes(self) -> int:
        return self.fea_solver.meshLoader.bdfInfo.nnodes

    @property
    def node_map(self) -> List[int]:
        """
        returns node map where you input a bdf node index and outputs the tacs node index
        """
        bdfNodes = range(self.num_nodes)
        return self.fea_solver.meshLoader.getLocalNodeIDsFromGlobal(bdfNodes,nastranOrdering=False)

    @property
    def function_names(self) -> List[str]:
        return self._function_names

    @property
    def analysis_dir(self) -> str:
        return self._analysis_dir

    @analysis_dir.setter
    def analysis_dir(self, new_dir):
        self._analysis_dir = new_dir

    @abstractmethod
    def __call__(self, dat_file:str, write_solution:bool=False):
        """
        the call method is required in all analysis functions to return funcs, sens
        """
        funcs = {}
        sens = {}
        return funcs, sens

    def f5_to_vtk(self):
        """
        get f5tovtk command
        """
        conversion_script = "~/git/tacs/extern/f5tovtk/f5tovtk"
        command = f"{conversion_script} {self.f5_filepath}"
        print(command)
        os.system(command)
        self._write_idx += 1
        


class MassStress(PytacsFunction):
    """
    Mass and stress Pytacs analysis function
    """
    def __init__(self, safety_factor:float=1.5, ks_weight:float=1000.0):
        super(MassStress,self).__init__(name="mass_stress")
        self._safety_factor = safety_factor
        self._ks_weight = ks_weight
        self._function_names = ['mass', 'ks_vmfailure']

        # for debugging
        self._nodes = None

    @property
    def nodes(self):
        assert(self._nodes is not None)
        return self._nodes

    def __call__(self, dat_file:str, write_solution:bool=False):

        #initialize pytacs with that data file
        self.fea_solver = pyTACS(dat_file)
            
        # Set up TACS Assembler
        self.fea_solver.initialize()

        #read the bdf & dat file into pytacs FEAsolver
        #SPs represents "StructuralProblems"
        SPs = self.fea_solver.createTACSProbsFromBDF()

        # Read in forces from BDF and create tacs struct problems
        for caseID in SPs:
            SPs[caseID].addFunction('mass', functions.StructuralMass)
            SPs[caseID].addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5, ksWeight=1000.0)

        func = {}; sens = {}
        for caseID in SPs:
            SPs[caseID].solve()
            SPs[caseID].evalFunctions(func,evalFuncs=self.function_names)
            SPs[caseID].evalFunctionsSens(sens,evalFuncs=self.function_names)
            self._nodes = SPs[caseID].getNodes()
            if write_solution:
                SPs[caseID].writeSolution(baseName=self.f5_filename, outputDir=self.analysis_dir)
                self.f5_to_vtk()
        
        return func, sens