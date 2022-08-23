
__all__ = ["CapsTacs"]

from capsManager.tacs import TacsAim
from capsManager.egads import EgadsAim
from caps2tacs.pytacs_functions import PytacsFunction
from typing import TYPE_CHECKING, List
from capsManager.tacs.design_variable import ThicknessVariable
import os


class CapsTacs:
    """
    Module to handle caps and tacs interface problems
    """
    def __init__(self, tacs_aim : TacsAim, egads_aim : EgadsAim, pytacs_function:PytacsFunction, compute_gradients:bool=True):
        """
        Module to handle caps and tacs interface problems
            tacs_aim : provide a fully setup tacs aim wrapper object
            egads_aim : provide a fully setup egads aim wrapper object
            pytacs_function : takes in dat_file and returns func, sens from TACS analysis
        """

        # make sure the aim wrapper modules are setup properly before proceeding
        assert(tacs_aim.is_setup)
        assert(egads_aim.is_setup)

        assert(isinstance(tacs_aim, TacsAim))
        assert(isinstance(egads_aim, EgadsAim))
        assert(isinstance(pytacs_function, PytacsFunction) or pytacs_function is None)

        # get the aim wrappers and pytacs function wrapper
        self._tacs_aim = tacs_aim
        self._egads_aim = egads_aim
        self._pytacs_function = pytacs_function

        # set the analysis directory of the pytacs function
        self.pytacs_function.analysis_dir = self.tacs_aim.analysis_dir

        # boolean of whether not to use derivatives
        self._compute_gradients = compute_gradients

        # func, gradient analysis counter to alternate analyses
        self._odd_analysis = True
        self._funcs = None
        self._sens = None

        # link the meshes
        self.tacs_aim.aim.input["Mesh"].link(self.egads_aim.aim.output["Surface_Mesh"])

    def sort_names(self, name_output:bool=True):
        name_list = [dv.name for dv in self.tacs_aim.design_variables]
        bool_list = [isinstance(dv,ThicknessVariable) for dv in self.tacs_aim.design_variables]
        zipped_list = zip(name_list, bool_list)
        sorted_ziplist = sorted(zipped_list, key = lambda x : x[0])
        if name_output:
            return [item[0] for item in sorted_ziplist]
        else: # otherwise bool output
            return [item[1] for item in sorted_ziplist]

    @property
    def design_variable_names(self) -> List[str]:
        return self.sort_names(name_output=True)

    @property
    def num_design_variables(self) -> int:
        return len(self.tacs_aim.design_variables)

    @property
    def list_is_struct_design_variable(self) -> List[bool]:
        return self.sort_names(name_output=False)

    @property
    def tacs_aim(self) -> TacsAim:
        assert(self._tacs_aim is not None)
        return self._tacs_aim

    @property
    def egads_aim(self) -> EgadsAim:
        assert(self._egads_aim is not None)
        return self._egads_aim

    @property
    def pytacs_function(self):
        assert(self._pytacs_function is not None)
        return self._pytacs_function

    @property
    def functions_dict(self) -> dict:
        assert(self._funcs is not None)
        return self._funcs

    @property
    def sensitivity_dict(self) -> dict:
        assert(self._sens is not None)
        return self._sens

    @property
    def function_names(self) -> List[str]:
        return self.pytacs_function.function_names

    @property
    def load_set_function_names(self) -> List[str]:
        return self.sensitivity_dict.keys()

    def set_function_name(self, function_name:str) -> str:
        found_match = False
        for name in self.load_set_function_names:
            found_match = function_name in name
            if found_match:
                return name
        if not found_match:
            return None

    def xpts_sensitivity(self, function_name:str):
        return self.sensitivity_dict[function_name]['Xpts']

    def struct_sensitivity(self, function_name:str):
        return self.sensitivity_dict[function_name]['struct']

    def reduced_analysis(self):
        """
        only compute the analysis every other time for function, gradient calls
        """
        if self._odd_analysis:
            self.analysis()
            self._odd_analysis = False
        else:
            self._odd_analysis = True

    def analysis(self, design_variables:List[float]=None, write_solution:bool=False):
        assert(self.pytacs_function is not None)
        
        if design_variables is not None:
            design_dict = self.design_variable_dictionary(design_variables)
            self.tacs_aim.update_design(design_dict)
            
        if not self.tacs_aim.is_setup:
            # if the tacs aim has been reconfigured then you have to setup the aim again
            self.tacs_aim.setup_aim()

        # run the pre analysis to generate dat file
        self.tacs_aim.pre_analysis()

        # run the pytacs function
        funcs, sens = self.pytacs_function(dat_file=self.tacs_aim.dat_file_path, write_solution=write_solution)

        # update funcs, sens dicts
        self._funcs = funcs

        # compute sensitivity product if requesting derivatives
        if self._compute_gradients:
            self._sens = sens
            self.coordinate_mesh_sensitivity_product()

        # force you to setup the tacs aim again next time
        self.tacs_aim._setup = False

    def design_variable_dictionary(self, design_variables:List[float]=None) -> dict:
        if design_variables is None:
            return None
        else:
            return {self.design_variable_names[idx] : design_variables[idx] for idx in range(self.num_design_variables)}

    def function(self, function_name:str, design_variables:List[float]=None):
        """
        return function values of TACS analysis
        """
        assert(function_name in self.function_names or function_name in self.load_set_function_names)
        self.analysis(design_variables)
        set_function_name = self.set_function_name(function_name=function_name)
        return self.functions_dict[set_function_name]

    def gradient(self, function_name:str, design_variables:List[float]=None):
        """
        Compute and report gradient values of TACS analysis for given function name
        returns grad_list, grad_dict
        """
        assert(function_name in self.function_names or function_name in self.load_set_function_names)
        self.analysis(design_variables)
        grad_dict = {}
        grad_list = []
        for dv_name in self.design_variable_names:
            set_function_name = self.set_function_name(function_name)
            grad_dict[dv_name] = self.tacs_aim.aim.dynout[set_function_name].deriv(dv_name)
            grad_list.append(grad_dict[dv_name])
        return grad_list, grad_dict

    def write_coordinates_test(self):
        """
        write the tacsAim sensitivity file to tacsAim analysis directory
        used for debugging and comparing the node indices of nastran and tacs for the  coordinate derivatives

        format of the .sens file:
            nnodes
            {for each node j}
            nodej xj yj zj
        """
        #where to print .sens file
        coordinates_file = os.path.join(self.tacs_aim.analysis_dir, "coordinate.out")

        nodal_coordinates = self.pytacs_function.nodes

        # node map from bdf to tacs nodes
        node_map = self.pytacs_function.node_map
        
        #open the file
        with open(coordinates_file, "w") as hdl:
            
            #for each node, print nodeind, dfdx, dfdy, dfdz for that mesh element
            for node_ind in range(self.pytacs_function.num_nodes): # d(Func1)/d(xyz)
                tacs_ind = node_map[node_ind]
                hdl.write("{} {} {} {}\n".format(node_ind + 1, nodal_coordinates[3*tacs_ind], nodal_coordinates[3*tacs_ind+1], nodal_coordinates[3*tacs_ind+2]))

    def write_sensitivity_file(self):
        """
        write the tacsAim sensitivity file to tacsAim analysis directory

        format of the .sens file:
            nfunc nDVsets

            {for each function i}
            func_i_key
            func_i_value
            nnodes
            nodej dfidxj dfidyj dfidzj
            {for each DVset k redundant for each function}
            DVset_k_key
            nDV_set_k
            dfunc_i_dx1 dfunc_i_dx2 ... dfunc_i_dxN (where N=nDV_set_k)
        """
        #where to print .sens file
        sensitivity_file = self.tacs_aim.sens_file_path

        # node map from bdf to tacs nodes
        node_map = self.pytacs_function.node_map
        
        #open the file
        with open(sensitivity_file, "w") as hdl:
            
            #write (nfunctions) in first line
            num_struct_DVs = sum(self.list_is_struct_design_variable)
            hdl.write(f"{self.pytacs_function.num_functions} {num_struct_DVs}\n")
            
            #for each function mass, stress, etc.
            for function_name in self.function_names:
                
                set_function_name = self.set_function_name(function_name=function_name)

                #get the pytacs/tacs sensitivity w.r.t. mesh for that function
                xpts_sens = self.xpts_sensitivity(set_function_name)
                
                #write the key,value,nnodes of the function
                hdl.write(f"{set_function_name}\n")
                hdl.write(f"{self.functions_dict[set_function_name]}\n")
                hdl.write(f"{self.pytacs_function.num_nodes}\n")

                #for each node, print nodeind, dfdx, dfdy, dfdz for that mesh element
                for node_ind in range(self.pytacs_function.num_nodes): # d(Func1)/d(xyz)
                    tacs_ind = node_map[node_ind]
                    hdl.write("{} {} {} {}\n".format(node_ind + 1, xpts_sens[3*tacs_ind], xpts_sens[3*tacs_ind+1], xpts_sens[3*tacs_ind+2]))

                # write the struct design variable derivatives, write each design variable separately
                struct_sens = self.struct_sensitivity(set_function_name)
                print(struct_sens)
                print(self.design_variable_names)
                print(self.list_is_struct_design_variable)
                for idx,is_struct_dv in enumerate(self.list_is_struct_design_variable):
                    if is_struct_dv:
                        struct_dv_name = self.design_variable_names[idx]
                        hdl.write(f"{struct_dv_name}\n")
                        hdl.write(f"{1}\n")
                        hdl.write(f"{struct_sens[idx]}\n")

    def coordinate_mesh_sensitivity_product(self):
        """
        write the tacsAim sensitivity file and run post analysis to multiply the 
        """
        self.write_sensitivity_file()
        self.tacs_aim.post_analysis()
        

        
