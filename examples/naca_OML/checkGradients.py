from __future__ import print_function

#from os import environ
#environ['CMPLX_MODE'] = "1"
from pyfuntofem.model  import *
from pyfuntofem.driver import *
from pyfuntofem.fun3d_interface import *

from pyoptsparse import SLSQP, Optimization
from mpi4py import MPI
import numpy as np

from pprint import pprint

from tacs import TACS, functions, constitutive, elements, pyTACS, problems
from pyfuntofem.tacs_interface import TacsSteadyInterface

#============================ Settings =======================================================#
#if complex off does adjoint, if complex on does complex step
#goal is to compare derivatives with adjoint and complex step

complexMode = True
if (complexMode): environ['CMPLX_MODE'] = "1"

analysis_type='aerothermoelastic'

#==============================TACS model setup==============================================#

class wedgeTACS(TacsSteadyInterface):
    def __init__(self, comm, tacs_comm, model, n_tacs_procs):
        super(wedgeTACS,self).__init__(comm, tacs_comm, model)

        assembler = None
        self.tacs_proc = False
        if comm.Get_rank() < n_tacs_procs:
            self.tacs_proc = True

            # Instantiate FEASolver
            structOptions = {
                'printtiming':True,
            }

            bdfFile = os.path.join(os.path.dirname(__file__), 'nastran_CAPS.dat')
            FEASolver = pyTACS(bdfFile, options=structOptions, comm=tacs_comm)

            # Material properties
            rho = 2780.0        # density kg/m^3
            E = 73.1e9          # Young's modulus (Pa)
            nu = 0.33           # Poisson's ratio
            kcorr = 5.0/6.0     # shear correction factor
            ys = 324.0e6        # yield stress
            specific_heat = 920.096
            cte = 24.0e-6
            kappa = 230.0

            tInput = 0.001*np.ones(3)

            def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
                elemIndex = kwargs['propID'] - 1
                t = tInput[elemIndex]

                prop = constitutive.MaterialProperties(rho=rho, specific_heat=specific_heat,
                                                       E=E, nu=nu, ys=ys, cte=cte, kappa=kappa)
                con = constitutive.IsoShellConstitutive(prop, t=t, tNum=dvNum)

                elemList = []
                transform = None
                for elemDescript in elemDescripts:
                    if elemDescript in ['CQUAD4', 'CQUADR']:
                        elem = elements.Quad4ThermalShell(transform, con)
                    else:
                        print("Uh oh, '%s' not recognized" % (elemDescript))
                    elemList.append(elem)

                # Add scale for thickness dv
                scale = [1.0]
                return elemList, scale

            # Set up elements and TACS assembler
            FEASolver.initialize(elemCallBack)
            assembler = FEASolver.assembler

        self._initialize_variables(assembler, thermal_index=6)
        self.initialize(model.scenarios[0],model.bodies)

    def post_export_f5(self):
        flag = (TACS.OUTPUT_CONNECTIVITY |
                TACS.OUTPUT_NODES |
                TACS.OUTPUT_DISPLACEMENTS |
                TACS.OUTPUT_STRAINS |
                TACS.OUTPUT_STRESSES |
                TACS.OUTPUT_EXTRAS)
        f5 = TACS.ToFH5(self.assembler, TACS.BEAM_OR_SHELL_ELEMENT, flag)
        f5.writeToFile('TACSoutput.f5')

#==================================================================================================#

maximum_mass = 40.0 
num_tacs_dvs = 3

# Set up the communicators
n_tacs_procs = 1

comm = MPI.COMM_WORLD
comm = comm

world_rank = comm.Get_rank()
if world_rank < n_tacs_procs:
    color = 55
    key = world_rank
else:
    color = MPI.UNDEFINED
    key = world_rank
tacs_comm = comm.Split(color,key)

#==================================================================================================#
# Originally _build_model()

thickness = 0.001

# Build the model
model = FUNtoFEMmodel('NACA Wing Simulation')
wing = Body('wing', analysis_type=analysis_type, group=0,boundary=1)

for i in range(num_tacs_dvs):
    wing.add_variable('structural',Variable('thickness '+ str(i),value=thickness,lower = 0.0001, upper = 0.01))

model.add_body(wing)

steady = Scenario('steady', group=0, steps=5)
function1 = Function('ksfailure',analysis_type='structural')
steady.add_function(function1)

function2 = Function('mass',analysis_type='structural',adjoint=False)
steady.add_function(function2)

model.add_scenario(steady)

#==================================================================================================#

# instantiate TACS on the master
solvers = {}
solvers['flow'] = Fun3dInterface(comm,model,flow_dt=1.0)
solvers['structural'] = wedgeTACS(comm,tacs_comm,model,n_tacs_procs)

# L&D transfer options
transfer_options = {'analysis_type': analysis_type,
                    'scheme': 'meld', 'thermal_scheme': 'meld'}

# instantiate the driver
driver = FUNtoFEMnlbgs(solvers,comm,tacs_comm,0,comm,0,transfer_options,model=model)
struct_tacs = solvers['structural'].assembler


#==================================================================================================#

# Comparison of derivatives with adjoint and complex step methods

tInput = 0.001*np.ones(3)
model.set_variables(tInput)

if (complexMode):

    #compute derivatives with complex step method

    x = []
    variables = self.model.get_variables()
    for i, var in enumerate(variables):
        if i == 0:
            x.append(var.value + 1j*epsilon)
        else:
            x.append(var.value)

    x = np.array(x)
    self.model.set_variables(x)

    self.driver.solve_forward()
    functions = self.model.get_functions()
    
    for index, func in enumerate(functions):
        print('Function %d Complex step: '%(index), func.value.imag/epsilon)

else:

    #compute derivatives with adjoint method
    
    self.driver.solve_forward()
    functions = self.model.get_functions()
    
    for index, func in enumerate(functions):
        print('Function %d Value'%(index), func.value)

    self.driver.solve_adjoint()

    grads = self.model.get_function_gradients()
    
    variables = self.model.get_variables()

    for i, func in enumerate(functions):
        for j, var in enumerate(variables):
            print("Adjoint Grad ", func.name, "Var: ", var.name, " ", grads[i][j])