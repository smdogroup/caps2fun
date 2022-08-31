
__all__ = ["PytacsFunction", "MassStress", "MassStressTransient"]

from typing import TYPE_CHECKING, List
from tacs import functions, pyTACS
from abc import ABC, abstractmethod
import os
import numpy as np

class PytacsFunction(ABC):
    """
    Base class for pytacs analysis functions
        can create your own custom analysis functions as well and use them in caps2tacs
        however need it to be subclass so that it will provide these properties -> coordinate derivatives, etc.
    """
    def __init__(self, name:str, transient:bool=False):
        self._name = name
        self._fea_solver = None
        self._function_names = []
        self._analysis_dir = None
        self._write_idx = 0
        self._load_set = 0

        self._transient = transient

    @property
    def load_set_str(self):
        set_str = str(self._load_set)
        rem_zeros = 3 - len(set_str)
        zero_str_list = ["0"] * rem_zeros
        zero_str = "".join(zero_str_list)
        return f"{zero_str}{set_str}"

    @property
    def name(self) -> str:
        return self._name

    @property
    def f5_base_filename(self) -> str:
        """
        name of f5 file without .f5 extension
        """
        return f"{self.name}_{self._write_idx}"

    @property
    def paraview_group_name(self) -> str:
        if self._transient:
            return f"{self.name}_{0}_{self.load_set_str}_..vtk"
        else:
            return f"{self.name}_{self.load_set_str}_..vtk"

    @property
    def f5_set_filename(self) -> str:
        return f"{self.f5_base_filename}_{self.load_set_str}"

    @property
    def f5_fixed_filename(self) -> str:
        return f"{self.name}_{self.load_set_str}_{self._write_idx}"

    @property
    def f5_fixed_filepath(self) -> str:
        return os.path.join(self.analysis_dir, f"{self.f5_fixed_filename}.f5")

    @property
    def f5_set_filepath(self) -> str:
        return os.path.join(self.analysis_dir, f"{self.f5_set_filename}.f5")

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
        get f5tovtk command done for each file
        ref for animation of a set of similar named paraview vtk files
        https://docs.paraview.org/en/latest/UsersGuide/dataIngestion.html
        """

        # first rename the f5 files to eliminate set extension
        conversion_script = "~/git/tacs/extern/f5tovtk/f5tovtk"
        if not self._transient:
            os.system(f"mv {self.f5_set_filepath} {self.f5_fixed_filepath}")
            command = f"{conversion_script} {self.f5_fixed_filepath}"
            print(command)
            os.system(command)
            self._write_idx += 1
        if self._transient:
            for filename in os.listdir(self.analysis_dir):
                if self.f5_base_filename in filename:
                    os.system(f"{conversion_script} {self.analysis_dir}/{filename}")
        


class MassStress(PytacsFunction):
    """
    Mass and stress Pytacs analysis function
        available function_names = [mass, ks_vmfailure]
    """
    def __init__(self, safety_factor:float=1.5, ks_weight:float=1000.0, steady:bool=True):
        super(MassStress,self).__init__(name="mass_stress", transient=False)
        self._safety_factor = safety_factor
        self._ks_weight = ks_weight
        self._function_names = ['mass', 'ks_vmfailure']
        self._steady = steady

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
            self._load_set = caseID-1
            SPs[caseID].solve()
            SPs[caseID].evalFunctions(func,evalFuncs=self.function_names)
            SPs[caseID].evalFunctionsSens(sens,evalFuncs=self.function_names)
            self._nodes = SPs[caseID].getNodes()
            if write_solution:
                SPs[caseID].writeSolution(baseName=self.f5_base_filename, outputDir=self.analysis_dir)
                self.f5_to_vtk()
        
        return func, sens

class MassStressTransient(PytacsFunction):
    """
    Mass and stress Pytacs analysis function
        available function_names = [mass, ks_vmfailure]
    """
    def __init__(self, safety_factor:float=1.5, ks_weight:float=1000.0, t0:float=0.0, tf:float=1.0, num_steps:int=100, amplitude=None):
        super(MassStressTransient,self).__init__(name="mass_stress", transient=True)
        self._safety_factor = safety_factor
        self._ks_weight = ks_weight
        self._t0 = t0
        self._tf = tf
        self._num_steps = num_steps
        self._amplitude = amplitude

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
        #SPs = self.fea_solver.createTACSProbsFromBDF()

        transientOptions = {'printlevel':1}
        TPs = self.fea_solver.createBDFtransientProb(tInit=self._t0, tFinal=self._tf, numSteps=self._num_steps, options=transientOptions, amplitude=self._amplitude)

        # Read in forces from BDF and create tacs struct problems
        for caseID in TPs:
            TPs[caseID].addFunction('mass', functions.StructuralMass)
            TPs[caseID].addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5, ksWeight=1000.0)

        func = {}; sens = {}
        for caseID in TPs:
            self._load_set = 0
            TPs[caseID].solve()
            TPs[caseID].evalFunctions(func,evalFuncs=self.function_names)
            TPs[caseID].evalFunctionsSens(sens,evalFuncs=self.function_names)
            self._nodes = TPs[caseID].getNodes()
            if write_solution:
                TPs[caseID].writeSolution(baseName=self.f5_base_filename, outputDir=self.analysis_dir)
                self.f5_to_vtk()
        
        return func, sens