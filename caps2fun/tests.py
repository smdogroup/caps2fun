#derivative testing, forward analysis, etc. for Caps2fun
import os
import numpy as np
from caps2fun import writeInput, readOutput

class Test():
    def __init__(self, DVdict, n_procs=192, functions=None):

        #copy the DV dict
        self.DVdict = DVdict

        #set the num procs
        self.n_procs = n_procs

        #set the functions to check
        self.functions = functions
        if (functions is None): 
            self.functions = ["ksfailure","cl","cd","mass"]

        #count the number of active DV
        self.nDV = 0
        for DV in DVdict:
            if (DV["active"]): self.nDV += 1

    def derivativeTest(self):

        #run the adjoint
        self.runAdjoint()

        #run the complex step
        self.runComplexStep()

        #compare the directional derivatives and write to file
        self.writeResults()

    def callCaps2fun(self):
        #call caps2fun through funtofem analysis

        #run funtofem analysis
        callMessage = "mpiexec_mpt -n {} python $CAPS2FUN/caps2fun/caps2fun.py 2>&1 > ./funtofem/output.txt".format(self.n_procs)
        os.system(callMessage)

    def runForward(self):
        #write the F2F input file for adjoint mode
        writeInput(self.DVdict, self.functions, mode="forward")

        #turnoff complex mode, this prob doesn't work
        #os.system("export CMPLX_MODE=0")

        #run funtofem
        self.callCaps2fun()

    def runAdjoint(self):
        #run the adjoint mode

        #write the F2F input file for adjoint mode
        writeInput(self.DVdict, self.functions, mode="adjoint")

        #turnoff complex mode, this prob doesn't work
        #os.system("export CMPLX_MODE=0")

        #run funtofem
        self.callCaps2fun()

        #read the output file
        self.adjoint_funcs, self.adjoint_grads = readOutput(self.DVdict, mode="adjoint")

    def runComplexStep(self,h=1e-30):
        #run funtofem in complex step mode

        #generate random perturbation for complex_step check
        x_dir = np.random.rand(self.nDV)
        x_dir = x_dir / np.linalg.norm(x_dir)

        #write the F2F input file for complex mode
        writeInput(self.DVdict, self.functions, mode="complex_step", eps=h, x_direction=x_dir)

        #setup complex mode, this prob doesn't work
        os.system("export CMPLX_MODE=1")

        #run funtofem
        self.callCaps2fun()

        #read the output file in complex mode, complex_step.out
        self.complex_funcs = readOutput(self.DVdict,mode="complex_step")

    def writeResults(self):
        #compare directional derivatives of adjoint vs complex step

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
        for function in self.functions:
            
            #write the function name "func,name"
            deriv_handle.write("func,{}\n".format(function))

            #write the adjoint and complex_step functions (real part) to the file
            deriv_handle.write("\tadjoint func      = {}\n".format(self.adjoint_funcs[function]))
            deriv_handle.write("\tcomplex_step func = {}\n".format(self.complex_funcs[function].real))

            #initialize directional derivative, and get this adjoint_grad
            adjoint_dderiv[function] = 0
            adjoint_grad = self.adjoint_grads[function]

            #loop over each design variable and dot product the perturbation and gradient to get directional derivative
            for i in range(self.nDV):
                adjoint_dderiv[function] += x_dir[i] * adjoint_grad[i]

            #write the adjoint directional derivative to the file
            line = "\tadjoint dderiv     = {}\n".format(adjoint_dderiv[function])
            deriv_handle.write(line)

            #compute complex step derivative, which is also for that direction
            #f(x+hp*1j)/h is the complex step
            complex_dderiv[function] = self.complex_funcs[function].imag

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