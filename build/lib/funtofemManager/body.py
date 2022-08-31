
__all__ = []

from typing import TYPE_CHECKING, List
from funtofemManager.analysis import Analysis
from funtofemManager.derivative import Derivative, Gradient
from funtofemManager.function import Function
from funtofemManager.variable import Variable, AerodynamicVariable, StructuralVariable

import numpy as np

MOTION_TYPES = ['deform', 'rigid', 'deform+rigid', 'rigid+deform']

class Body:
    """
    FUNtoFEM body class, is used to represent deformable bodies or aircraft structures in the flow.
    """
    def __init__(self, name:str, analysis:Analysis, id:int=0, group:int=0,
                boundary:int=0, use_fun3d:bool=True, motion_type:str='deform'):
        """
        Parameters
        ----------
        name - name of body
        id - ID number of the body in the model/problem
        group - group # of the body used for coupled bodies
        boundary - the FUN3D boundary number associated with the body
        use_fun3d - whether to use FUN3D for the analysis
        motion_type - 
        """
        assert(motion_type in MOTION_TYPES)

        self._name = name
        self._analysis = analysis
        self._id = id
        self._group = group
        self._boundary = boundary
        self._use_fun3d = use_fun3d
        self._motion_type = motion_type

        # variables and derivatives
        self._variables_are_setup = False
        self._variables = []
        self._gradients = []

        # transfer booleans
        self._elastic_transfer = None
        self._thermal_transfer = None

        # number of nodes
        self._struct_nnodes = None
        self._aero_nnodes = None

        # number of degrees of freedom
        self._elastic_ndof = 3
        self._thermal_ndof = 1

        # thermal settings
        self._thermal_index = 3
        self._temp_ref = 300.0 # Kelvin

        # node ids
        self._struct_ids = None
        self._aero_ids = None

        # structural elastic transfer variables
        self._struct_disps = None
        self._struct_forces = None
        self._struct_loads = None
        self._struct_shape_term = None

        # aerodynamic elastic transfer variables
        self._rigid_transform = None
        self._aero_disps = None
        self._aero_forces = None
        self._aero_loads = None
        self._aero_shape_term = None

        # structural thermal transfer variables
        self._struct_temps = None
        self._struct_heat_flux = None

        # aerodynamic thermal transfer variables
        self._aero_temps = None
        self._aero_heat_flux = None

    @property
    def id(self) -> int:
        return self._id

    @id.setter
    def id(self, new_id:int):
        self._id = new_id

    @property
    def variables(self) -> List[Variable]:
        return self._variables

    @property
    def aerodynamic_variables(self) -> List[AerodynamicVariable]:
        return [var for var in self.variables if isinstance(var, AerodynamicVariable)]

    @property
    def structural_variables(self) -> List[StructuralVariable]:
        return [var for var in self.variables if isinstance(var, StructuralVariable)]

    @property
    def num_active_variables(self) -> int:
        return len(self.active_variables)

    @property
    def active_variables(self) -> List[Variable]:
        return [var for var in self.variables if var.active]

    def add_variable(self, new_variable:Variable):
        self._variables += [new_variable]

    def add_function(self, new_function:Function):
        """
        adds a function tracked by the body and thus will add a gradient for this
        """

        # ensures all the variables are setup before adding functions
        assert(self.variables_are_setup)

        # make a list of derivatives for the new 
        derivative_list = []
        for variable in self.variables:
            derivative_list.append(Derivative(function=new_function, variable=variable))

        # setup a gradient object for this body and this function and add to list of gradients
        new_gradient = Gradient(function=function, derivatives=derivative_list)
        self._gradients.append(new_gradient)

    @property
    def variables_are_setup(self) -> bool:
        return self._variables_are_setup

    @variables_are_setup.setter
    def variables_are_setup(self, new_bool:bool):
        self._variables_are_setup = new_bool

    def collect_coordinate_derivatives(self, comm, discipline:str, root:int=0):
        """
        Used for writing the sensitivity files for the aerodynamic and structural 
        meshes on the root processor. The code collects sensitivities from each proc
        and collects the results on the root node
        """

        assert(discipline in ["aerodynamic", "structural"])

        if discipline == "aerodynamic":
            all_aero_ids = comm.gather(self._aero_ids, root=root)
            all_aero_shape = comm.gather(self._aero_shape_term, root=root)

            aero_ids = []
            aero_shape = []

            if comm.rank == root:
                # Discard any entries that are None
                aero_ids = []
                for d in all_aero_ids:
                    if d is not None:
                        aero_ids.append(d)

                aero_shape = []
                for d in all_aero_shape:
                    if d is not None:
                        aero_shape.append(d)

                if len(aero_shape) > 0:
                    aero_shape = np.concatenate(aero_shape)
                else:
                    aero_shape = np.zeros((3, 1))

                if len(aero_ids) == 0:
                    aero_ids = np.arange(aero_shape.shape[0]//3, dtype=int)
                else:
                    aero_ids = np.concatenate(aero_ids)

            return aero_ids, aero_shape

        elif discipline == "structural":
            all_struct_ids = comm.gather(self._struct_ids, root=root)
            all_struct_shape = comm.gather(self._struct_shape_term, root=root)

            struct_ids = []
            struct_shape = []

            if comm.rank == root:
                # Discard any entries that are None
                struct_ids = []
                for d in all_struct_ids:
                    if d is not None:
                        struct_ids.append(d)

                struct_shape = []
                for d in all_struct_shape:
                    if d is not None:
                        struct_shape.append(d)

                if len(struct_shape) > 0:
                    struct_shape = np.concatenate(struct_shape)
                else:
                    struct_shape = np.zeros((3, 1))

                if len(struct_ids) == 0:
                    struct_ids = np.arange(struct_shape.shape[0]//3, dtype=int)
                else:
                    struct_ids = np.concatenate(struct_ids)

            return struct_ids, struct_shape
