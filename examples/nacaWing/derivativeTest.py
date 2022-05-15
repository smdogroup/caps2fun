#derivative testing for runF2F
import os
import numpy as np
from fileIO import writeInput, makeDVdict, readOutput

#set the num procs
n_procs = 192

#analysis type, "aerothermal", "aeroelastic", "aerothermoelastic"
f2fanalysis = "aeroelastic"

#fun3d analysis, "inviscid", "laminar", "turbulent"
fun3danalysis = "inviscid"

#set the functions to check
functions = ["ksfailure","cl","cd","mass"]

#----------------------------------------------------------------------#

#generate variables and functions from writeInput.py
DVdict = makeDVdict()

#count the number of active DV
nDV = 0
for DV in DVdict:
    if (DV["active"]): nDV += 1

print("Initialized DV dicts\n", flush=True)

#----------------   Run FUNtoFEM in Complex mode ----------------------#

#generate random perturbation for complex_step check
h = 1.0e-30
x_dir = np.random.rand(nDV)
x_dir = x_dir / np.linalg.norm(x_dir)
print("Generated perturbation\n", flush=True)

#write the F2F input file for complex mode
print("Trying to write input...", flush=True)
writeInput(DVdict, functions, f2fanalysis=f2fanalysis, fun3danalysis=fun3danalysis, mode="complex_step", eps=h, x_direction=x_dir)
print("Wrote input\n", flush=True)

#setup complex mode, this prob doesn't work
os.system("export CMPLX_MODE=1")

#run funtofem
print("Trying to run F2F in complex mode...", flush=True)
callMessage2 = "mpiexec_mpt -n {} python runF2F.py 2>&1 > ./funtofem/output.txt".format(n_procs)
os.system(callMessage2)
print("Ran F2F in complex mode\n", flush=True)

#read the output file in complex mode, complex_step.out
complex_funcs = readOutput(DVdict,mode="complex_step")

#-----------------  Run FUNtoFEM with Adjoint mode --------------------#

#write the F2F input file for adjoint mode
writeInput(DVdict, functions, f2fanalysis=f2fanalysis, fun3danalysis=fun3danalysis,  mode="adjoint")

#turnoff complex mode, this prob doesn't work
os.system("export CMPLX_MODE=0")

#run funtofem
callMessage = "mpiexec_mpt -n {} python runF2F.py 2>&1 > ./funtofem/output.txt".format(n_procs)
os.system(callMessage)

#read the output file
adjoint_funcs, adjoint_grads = readOutput(DVdict, mode="adjoint")

#----------------  Compare Results of Adjoint and Complex Step -------------------------------------#

#make a derivative_check.out file in funtofem folder
funtofemFolder = os.path.join(os.getcwd(), "funtofem")
deriv_file = os.path.join(funtofemFolder, "derivative_check.out")

deriv_handle = open(deriv_file, "w")

fileExists = os.path.exists(deriv_file)

#write a an introductory line to the file
if (not(fileExists)):
    deriv_handle.write("Funtofem Derivative Check, Adjoint vs Complex Step\n")
    deriv_handle.write("----------------------------")
    deriv_handle.write("----------------------------\n")

adjoint_dderiv = {}
complex_dderiv = {}

#for each function in the analysis
for function in functions:
    
    #write the function name "func,name"
    deriv_handle.write("func,{}\n".format(function))

    #write the adjoint and complex_step functions (real part) to the file
    deriv_handle.write("\tadjoint func      = {}\n".format(adjoint_funcs[function]))
    deriv_handle.write("\tcomplex_step func = {}\n".format(complex_funcs[function].real))

    #initialize directional derivative, and get this adjoint_grad
    adjoint_dderiv[function] = 0
    adjoint_grad = adjoint_grads[function]

    #loop over each design variable and dot product the perturbation and gradient to get directional derivative
    for i in range(nDV):
        adjoint_dderiv[function] += x_dir[i] * adjoint_grad[i]

    #write the adjoint directional derivative to the file
    line = "\tadjoint dderiv     = {}\n".format(adjoint_dderiv[function])
    deriv_handle.write(line)

    #compute complex step derivative, which is also for that direction
    #f(x+hp*1j)/h is the complex step
    complex_dderiv[function] = complex_funcs[function].imag

    #write the complex step directional derivative to the file
    line = "\tcomplex_step dderiv = {}\n".format(complex_dderiv[function])
    deriv_handle.write(line)

    #compute the absolute error
    absoluteError = abs( complex_dderiv[function] - adjoint_dderiv[function] )

    #write the absolute error to the file
    line = "\tabsolute error      = {}\n".format(absoluteError)
    deriv_handle.write(line)

    #compute the relative error
    if (abs(complex_dderiv[function]) < 1.0e-15):
        relativeError = None
    else:
        relativeError = absoluteError / abs(complex_dderiv[function])

    #write relative error to a file
    line = "\relative error       = {}\n".format(relativeError)
    deriv_handle.write(line)


#close the derivative_check.out file
deriv_handle.close()