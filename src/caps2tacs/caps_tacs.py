
__all__ = ["CapsTacs"]

from capsManager.tacs import TacsAim
from capsManager.egads import EgadsAim
from tacsManager.pytacs_functions import PytacsFunction
from typing import TYPE_CHECKING, List
from capsManager.tacs.design_variable import ThicknessVariable
import os, sys
import matplotlib.pyplot as plt


class CapsTacs:
    """
    Module to handle caps and tacs interface problems
    """
    def __init__(self, name:str, comm, tacs_aim : TacsAim, egads_aim : EgadsAim, pytacs_function:PytacsFunction, 
        compute_gradients:bool=True, write_solution:bool=True, view_plots:bool=False, report_history:bool=False):
        """
        Module to handle caps and tacs interface problems
            tacs_aim : provide a fully setup tacs aim wrapper object
            egads_aim : provide a fully setup egads aim wrapper object
            pytacs_function : takes in dat_file and returns func, sens from TACS analysis
        """

        self._name = name
        self._comm=comm

        # make sure the aim wrapper modules are setup properly before proceeding
        if self.root_proc:
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
        if pytacs_function is not None:
            self.pytacs_function.analysis_dir = self.analysis_dir

        # boolean of whether not to use derivatives
        self._compute_gradients = compute_gradients
        self._write_solution = write_solution

        # func, gradient analysis counter to alternate analyses
        self._odd_analysis = True
        self._funcs = None
        self._sens = None

        # boolean to view paraview group on destructor
        self._view_plots = view_plots
        self._report_history = report_history
        self._function_history = {}

        # link the meshes
        if self.root_proc:
            self.tacs_aim.aim.input["Mesh"].link(self.egads_aim.aim.output["Surface_Mesh"])

    @property
    def root_proc(self) -> bool:
        return self._comm is None or self._comm.Get_rank() == 0

    @property
    def multi_proc(self) -> bool:
        return self._comm is not None

    @property
    def analysis_dir(self) -> str:
        the_dir = None
        if self.root_proc:
            the_dir = self.tacs_aim.analysis_dir
        if self.multi_proc:
            the_dir = self._comm.bcast(the_dir,root=0)
        return the_dir

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
    def name(self) -> str:
        return self._name

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

    def analysis_gatekeeper(self, design_dict:dict=None) -> bool:
        """
        gatekeeper function for the analysis, that ensures to not rerun the analysis when the design hasn't changed
        returns bool whether the design was changed or not
        """
        changed_design = False
        if design_dict is not None:
            #design_dict = self.design_variable_dictionary(design_variables)
            if self.root_proc:
                changed_design = self.tacs_aim.update_design(design_dict)
            if self.multi_proc:
                changed_design = self._comm.bcast(changed_design, root=0)
            if changed_design:
                #self.tacs_aim._geometry.view()
                self.analysis()
        
        else: # run the analysis again if not design dict
            if self._funcs is None or self._sens is None:
                self.analysis()
                changed_design = True # probably first design

        return changed_design

    def analysis(self):
        assert(self.pytacs_function is not None)

        if not self.tacs_aim.is_setup and self.root_proc:
            # if the tacs aim has been reconfigured then you have to setup the aim again
            self.tacs_aim.setup_aim()

        # run the pre analysis to generate dat file
        if self.root_proc:
            self.tacs_aim.pre_analysis()

        # run the pytacs function
        funcs, sens = self.pytacs_function(dat_file=self.tacs_aim.dat_file_path, write_solution=self._write_solution)

        # update funcs, sens dicts
        self._funcs = funcs

        # compute sensitivity product if requesting derivatives
        if self._compute_gradients:
            self._sens = sens
            self.coordinate_mesh_sensitivity_product()

        # force you to setup the tacs aim again next time
        if self.root_proc:
            self.tacs_aim._setup = False

    def design_variable_dictionary(self, design_variables:List[float]=None) -> dict:
        if design_variables is None:
            return None
        else:
            return {self.design_variable_names[idx] : design_variables[idx] for idx in range(self.num_design_variables)}

    def function(self, function_name:str, design_dict:dict=None):
        """
        return function values of TACS analysis
        """
        assert(function_name in self.function_names or function_name in self.load_set_function_names)
        self.analysis_gatekeeper(design_dict)

        set_function_name = self.set_function_name(function_name=function_name)
        function_value = self.functions_dict[set_function_name]

        self.register_history(name=function_name, value=function_value)
        return function_value

    def register_history(self, name:str, value:float):
        if not (name in self._function_history):
            self._function_history[name] = []
        self._function_history[name].append(value)

    def gradient(self, function_name:str, design_dict:dict=None):
        """
        Compute and report gradient values of TACS analysis for given function name
        returns grad_list, grad_dict
        """
        assert(function_name in self.function_names or function_name in self.load_set_function_names)
        self.analysis_gatekeeper(design_dict)
        grad_dict = {}
        #grad_list = []
        if self.root_proc:
            for dv_name in design_dict:
                load_set_function_name = self.set_function_name(function_name)
                grad_dict[dv_name] = self.tacs_aim.aim.dynout[load_set_function_name].deriv(dv_name)
                #grad_list.append(grad_dict[dv_name])
        if self.multi_proc:
            grad_dict = self._comm.bcast(grad_dict,root=0)
        return grad_dict 

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
                #print(struct_sens)
                #print(self.design_variable_names)
                #print(self.list_is_struct_design_variable)
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
        if self.root_proc:
            self.write_sensitivity_file()
            self.tacs_aim.post_analysis()

    def view_paraview_group(self):
        os.system(f"paraview {self.tacs_aim.analysis_dir}/{self.pytacs_function.paraview_group_name}")

    def plot_function_history(self):
        """
        plot the function history and save it during the __del__ step
        had this function in the __del__ destructor but doesn't work very well bc matplotlib doesn't work in this __del__ stage
        """
        history = self._function_history
        #print(history)
        num_iters = min([len(history[key]) for key in history])
        max_vals = {key:max(history[key]) for key in history}
        if num_iters > 0:
            #fig = plt.figure("CapsTacs")
            iterations = [idx+1 for idx in range(num_iters)]
            colors="kbg"
            color_ind = 0
            for key in history:
                values = history[key][:num_iters]
                norm_values = [value/max_vals[key] for value in values]
                style = colors[color_ind] + "-"
                color_ind += 1
                label=f"{key}/{max_vals[key]:.2f}"
                print(iterations, norm_values)
                plt.plot(iterations, norm_values, style, label=label, lw=3)
                
            plt.legend()
            plt.xlabel("iterations")
            plt.ylabel("functions")
            plt.title(f"Caps2tacs optimization of '{self.name}'")
            
            #plt.show()
            filepath = os.path.join(self.tacs_aim.analysis_dir, f"{self.name}_opt.png")
            plt.savefig(filepath)
        else:
            raise AssertionError("Can't plot function history when no iterations have been recorded yet...")

    def print_function_history(self):
        history = self._function_history
        #print(history)
        num_iters = min([len(history[key]) for key in history])
        if num_iters > 0:
            for key in history:
                values = history[key][:num_iters]
                print(f"{key} history:")
                print(f"\t{values}")
        else:
            raise AssertionError("Can't plot function history when no iterations have been recorded yet...")

    def print_final_design(self):
        """
        TODO : make it save the best design from each step
        for now just outputs the last design
        """
        optimal_design_dict = {dv.name : dv.value for dv in self.tacs_aim.design_variables}
        print(f"optimal design = {optimal_design_dict}")

    def __del__(self):
        """
        destructor for caps to tacs problem, tells how to view results after done
        """

        if self._report_history and self.root_proc:
            self.print_function_history()
            self.print_final_design()
            
        print(f"\nYou've now finished your caps2tacs analysis/optimization")
        
        if self._view_plots:
            print(f"\tYou can find .vtk files for paraview in the following directory and use 'paraview' to open the batch of files")
            if self.root_proc:
                print(f"\tcd {self.tacs_aim.analysis_dir}")
                print(f"\tparaview {self.pytacs_function.paraview_group_name}")
            print(f"Once inside paraview the following command uses the u,v,w displacement field to apply deformation on the animation batch")
            print(f"\tu*iHat+v*jHat+w*kHat")
            print(f"Once you save the deformation results or field output animations (have to save as group pngs) go to the following website to make a gif...")
            print(f"\thttps://www.freeconvert.com/png-to-gif")
            print(f"Then return to the orig directory with 'cd -'")

            self.view_paraview_group()

        

        

        
