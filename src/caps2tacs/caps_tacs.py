
__all__ = ["CapsTacs"]

from capsManager.tacs import TacsAim
from capsManager.egads import EgadsAim
from caps2tacs.pytacs_functions import PytacsFunction
from typing import TYPE_CHECKING, List


class CapsTacs:
    """
    Module to handle caps and tacs interface problems
    """
    def __init__(self, tacs_aim : TacsAim, egads_aim : EgadsAim, pytacs_function:PytacsFunction):
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

        # func, gradient analysis counter to alternate analyses
        self._odd_analysis = True
        self._funcs = None
        self._sens = None

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

    def xpts_sensitivity(self, function_name:str):
        return self.sensitivity_dict[function_name]['Xpts']

    def struct_sensitivity(self, function_name:str):
        return self.sensitivity_dict[function_name]['struct']

    def analysis(self, write_solution:bool=False):
        assert(self.pytacs_function is not None)

        if self._odd_analysis:
            # link the meshes
            self.tacs_aim.aim.input["Mesh"].link(self.egads_aim.aim.output["Surface_Mesh"])

            # run the pre analysis to generate dat file
            self.tacs_aim.pre_analysis()

            # run the pytacs function
            funcs, sens = self.pytacs_function(dat_file=self.tacs_aim.dat_file_path, write_solution=write_solution)

            # update funcs, sens dicts
            self._funcs = funcs
            self._sens = sens

            self._odd_analysis = False
        else:
            self._odd_analysis = True

    def function(self, function_name:str):
        assert(function_name in self.function_names)
        return self.functions_dict[function_name]

    def coordinate_mesh_sensitivity_product(self):
        #where to print .sens file
        sensitivity_file = self.tacs_aim.sens_file_path

        # node map from 
        node_map = self.pytacs_function.node_map
        
        #open the file
        with open(sensitivity_file, "w") as f:
            
            #write (nfunctions) in first line
            f.write("{}\n".format(self.pytacs_function.num_functions))
            
            #for each function mass, stress, etc.
            for function_name in self.function_names:
                
                #get the pytacs/tacs sensitivity w.r.t. mesh for that function
                xpts_sens = self.xpts_sensitivity(function_name)
                
                #write the key,value,nnodes of the function
                f.write(function_name + "\n")
                f.write("{}\n".format(self.func[function_name]))
                f.write("{}\n".format(self.pytacs_function.num_nodes))

                #for each node, print nodeind, dfdx, dfdy, dfdz for that mesh element
                for node_ind in range(self.pytacs_function.num_nodes): # d(Func1)/d(xyz)
                    bdf_ind = node_ind + 1
                    tacs_ind = node_map[bdf_ind]
                    f.write("{} {} {} {}\n".format(bdf_ind, xpts_sens[tacs_ind,0], xpts_sens[tacs_ind,1], xpts_sens[tacs_ind,2]))
        
    def compute_shape_derivatives(self):
        self.coordinate_mesh_sensitivity_product()
        self.tacs_aim.post_analysis()
        
        # store derivatives
        # for key in self.funcKeys:
            
        #     self.func[key] = self.tacsAim.dynout[key].value
        #     self.grad[key] = np.zeros((self.nvar))

        #     print("Function {} = {:5f}".format(key,self.func[key]))

        #     #loop over each design variable to get the full df/dD gradient
        #     #print("struct sens: ",self.sens[key]['struct'])
        #     ind = 0
        #     thickind = 0
        #     for desvar in self.desvars:
        #         if ("thick" in desvar):
        #             #use struct here, #print(self.sens[key]['struct'])
        #             #struct includes geomDVs in same order
        #             self.grad[key][ind] = self.sens[key]['struct'][ind]
        #             thickind += 1
        #         else:
        #             self.grad[key][ind] = self.tacsAim.dynout[key].deriv(desvar)

        #         print("\tdf/d{} = {:5f}".format(desvar, self.grad[key][ind]))
        #         ind += 1

    def gradient(self, function_name:str):
        assert(function_name in self.function_names)

        # compute the coordinate mesh sensitivity product
        

        
