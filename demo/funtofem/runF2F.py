#!/usr/bin/env python
"""
This file is part of the package FUNtoFEM for coupled aeroelastic simulation
and design optimization.

Copyright (C) 2015 Georgia Tech Research Corporation.
Additional copyright (C) 2015 Kevin Jacobson, Jan Kiviaho and Graeme Kennedy.
All rights reserved.

FUNtoFEM is licensed under the Apache License, Version 2.0 (the "License");
you may not use this software except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

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



#==================================================================================================#
# Originally tacs_model.py

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

            t = 0.001

            def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
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

analysis_type='aerothermoelastic'

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

obj_scale = 0.0106
con_scale = 3.35
    


def objFunc(self, xdict):

    tInput1 = xdict["xvars"]

    model.set_variables(tInput1)
    driver.solve_forward()
    functions = model.get_functions()

    func1 = functions[0].value * 1.0/obj_scale
    func2 = functions[1].value * 1.0/con_scale

    funcs = {}
    funcs["obj"] = func1
    funcs["con"] = func2

    fail = False

    if comm.rank == 0:
        print('\n---------- FUNCTION SOLVE ----------', flush=True)
        print('\nDesign Vars:       ', tInput1, flush=True) 
        print('\nObjective Value:   ', func1, flush=True)
        print('\nConstraint Value:  ', func2, flush=True)

        print('\n---------- FUNCTION SOLVE ----------', flush=True, file=optHist)
        print('\nDesign Vars:       ', tInput1, flush=True, file=optHist) 
        print('\nObjective Value:   ', func1, flush=True, file=optHist)
        print('\nConstraint Value:  ', func2, flush=True, file=optHist)

    print('\ncomm.rank:      ', comm.rank, flush=True, file=optHistAll)
    print('\n---------- FUNCTION SOLVE ----------', flush=True, file=optHistAll)
    print('\nDesign Vars:       ', tInput1, flush=True, file=optHistAll) 
    print('\nObjective Value:   ', func1, flush=True, file=optHistAll)
    print('\nConstraint Value:  ', func2, flush=True, file=optHistAll)

    return funcs, fail


def objGrad(self, xdict, funcs): 

    tInput1 = xdict["xvars"]

    model.set_variables(tInput1)
    driver.solve_adjoint()
    grads = model.get_function_gradients()

    grad1 = np.array(grads[0][:]) * 1.0/obj_scale
    grad2 = np.array(grads[1][:]) * 1.0/con_scale

    sens = {}
    sens = {
        "obj": {
            "xvars": [grad1]
        },
        "con": {
            "xvars": [grad2]
        },
    }

    fail = False

    if comm.rank == 0:
        print('\n---------- GRADIENT SOLVE ----------', flush=True)
        print('\nDesign Vars:            ', tInput1, flush=True) 
        print('\nObjective Grad Value:   ', grad1, flush=True)
        print('\nConstraint Grad Value:  ', grad2, flush=True)

        print('\n---------- GRADIENT SOLVE ----------', flush=True, file=optHist)
        print('\nDesign Vars:            ', tInput1, flush=True, file=optHist) 
        print('\nObjective Grad Value:   ', grad1, flush=True, file=optHist)
        print('\nConstraint Grad Value:  ', grad2, flush=True, file=optHist)
    
    print('\ncomm.rank:           ', comm.rank, flush=True, file=optHistAll)
    print('\n---------- GRADIENT SOLVE ----------', flush=True, file=optHistAll)
    print('\nDesign Vars:            ', tInput1, flush=True, file=optHistAll) 
    print('\nObjective Grad Value:   ', grad1, flush=True, file=optHistAll)
    print('\nConstraint Grad Value:  ', grad2, flush=True, file=optHistAll)

    return sens, fail





# dp = wedge_adjoint(analysis_type='aerothermoelastic') # 'aeroelastic') # 'aerothermoelastic') # 'aerothermal')

tInput = 0.001*np.ones(3)
model.set_variables(tInput)

driver.solve_forward()
functionEvals = model.get_functions()
print(functionEvals)

driver.solve_adjoint()
grads = model.get_function_gradients()
print(grads)



# optProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", dp.objFunc)

# optProb.addVarGroup("xvars", 12, "c", lower=0.0001*np.ones(12), upper=0.01*np.ones(12), value=0.001)

# optProb.addConGroup("con", 1, lower=1, upper=1)
# optProb.addObj("obj")

# comm = MPI.COMM_WORLD
# # if comm.rank == 0:
# print(optProb)

# optOptions = {"IPRINT": -1}
# opt = SLSQP(options=optOptions)
# sol = opt(optProb, sens=dp.objGrad)

# if comm.rank == 0:
#     print(sol)
#     print(sol, file=dp.optHist)
#     print('\nsol.xStar:  ', sol.xStar)
#     print('\nsol.xStar:  ', sol.xStar, file=dp.optHist)
# print(sol, file=dp.optHistAll)
# print('\nsol.xStar:  ', sol.xStar, file=dp.optHistAll)

# dp.optHist.close()
# dp.optHistAll.close()
