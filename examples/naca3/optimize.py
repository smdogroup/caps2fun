'''
script to run funtofem from design variable inputs, funtofem.input file and .csm files
returns function, gradient information in funtofem.out file
'''
from pyoptsparse import SLSQP, Optimization
import os
import numpy as np

class WingOptimizer():
    #class to run pyoptsparse optimize on the outside of funtofem

    def __init__(self, DVdict):
        #set the DV dict here
        self.DVdict = DVdict

        #other information such as analysis type, mesher type is  setup in funtofem.py
        self.n_procs = 192

        #make status file
        statusFile = os.path.join(os.getcwd(), "status.txt")
        self.status = open(statusFile, "w")

        #iterations
        self.iteration = 1

        self.cwrite("Aerothermoelastic Optimization with FuntoFem and ESP/CAPS\n")
        self.cwrite("\tDesign Problem: NACA Symmetric Wing\n")
        self.cwrite("Authors: Sean Engelstad, Sejal Sahu, Graeme Kennedy\n")
        self.cwrite("\tGeorgia Tech SMDO Lab April 2022\n")
        self.cwrite("----------------------------")
        self.cwrite("----------------------------\n")

    def cwrite(self, msg):
        self.status.write(msg)
        self.status.flush()

    def writeInput(self, x):
        # script to write the design variables to funtofem
        # sends them in the file funtofem.in
        # also tells which kind of analysis to perform, etc.

        #make sure funtofem folder exists
        funtofemFolder = os.path.join(os.getcwd(), "funtofem")
        if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

        #make funtofem input file
        inputFile = os.path.join(funtofemFolder, "funtofem.in")
        inputHandle =  open(inputFile, "w")

        #overwrite values in dict with the input x
        structDVs = x["struct"]
        shapeDVs = x["shape"]

        structCt = 0
        shapeCt = 0
        for DV in self.DVdict:
            #get attributes
            name = DV["name"]
            dvType = DV["type"]
            capsGroup = DV["capsGroup"]

            #get value
            if (dvType == "shape"):
                value = shapeDVs[shapeCt]
                shapeCt += 1
            elif (dvType == "struct"):
                value = structDVs[structCt]
                structCt += 1

            #make the line strings
            if (dvType == "shape"):
                line = "{},{},{}\n".format(name, dvType, value)
            elif (dvType == "struct"):
                line = "{},{},{},{}\n".format(name, dvType, capsGroup, value)   

            #write the line to the input file
            inputHandle.write(line)

        #close the input file
        inputHandle.close()

        #update status to the inputs that were run
        self.cwrite("Global Iteration #{}\n".format(self.iteration))
        self.cwrite("\tShape DVs {}\n".format(shapeDVs))
        self.cwrite("\tStruct DVs {}\n".format(structDVs)) 
        self.iteration += 1           

    def callFuntoFem(self):
        #call funtofem via a system call (os.system)
        #the reason for this is fun3d can't be run twice in the system python script
        #this gets around that issue

        #update status
        self.cwrite("\tRunning F2F... ")

        #instead of bash, do system call inside of this python script
        callMessage = "mpiexec_mpt -n {} python runF2F.py 2>&1 > output.txt".format(self.n_procs)
        os.system(callMessage)

        #update status that F2F finished
        self.cwrite("finished F2F!\n")

    def readOutput(self):
        #read functions, gradients from funtofem call
        #make sure funtofem folder exists
        funtofemFolder = os.path.join(os.getcwd(), "funtofem")
        if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

        #make funtofem input file
        outputFile = os.path.join(funtofemFolder, "funtofem.out")
        outputHandle =  open(outputFile, "r")

        lines = outputHandle.readlines()

        #read first line
        firstline = lines[0]
        parts = firstline.split(",")
        nfunc = int(parts[0])
        nshape = int(parts[1])
        nstruct = int(parts[2])

        self.functions = {}
        self.shapeGrad = {}
        self.structGrad = {}

        #read function values and gradients
        ifunc = -1
        for line in lines:
            if ("func" in line):
                #func, name, value
                ifunc += 1
                parts = line.split(",")
                name = parts[1]
                value = float(parts[2])

                self.functions[name] = value
                self.shapeGrad[name] = np.zeros((nshape))
                self.structGrad[name] = np.zeros((nstruct))
            
            elif ("grad" in line):
                #grad,* line
                if ("shape" in line):
                    #grad,shape,deriv1,deriv2,...,derivN
                    parts = line.split(",")
                    parts = parts[2:]
                    for ishape in range(len(parts)):
                        self.shapeGrad[name][ishape] = float(parts[ishape])
                elif ("struct" in line):
                    #grad,struct,deriv1,deriv2,...,derivN
                    parts = line.split(",")
                    parts = parts[2:]
                    for istruct in range(len(parts)):
                        self.structGrad[name][istruct] = float(parts[istruct])

        outputHandle.close()

        #update status
        self.cwrite("\tRead output files for func, sens\n")

    def objCon(self, x):
        #py opt sparse function evaluator

        #write input, call funtofem, and get results
        self.writeInput(x)
        self.callFuntoFem()
        self.readOutput()

        #objective, constraint functions
        funcs = {}
        funcs["obj"] = self.functions["ksfailure"]
        #funcs["con"] = self.functions["mass"]

        fail = False

        return funcs, fail

    def objGrad(self, x, funcs):
        
        #read output just in case
        self.readOutput()
        shapeGrad = self.shapeGrad["ksfailure"]
        structGrad = self.structGrad["ksfailure"]

        #assume funtofem and esp/caps grid generation already done in objCon call
        fail = False

        sens = {}
        sens["obj"] = { "struct": shapeGrad,
                        "shape" : massShstructGradapeGrad}
        #sens["con"] = {"struct" : 0 * massStructGrad,
        #                "shape" : 0 * massShapeGrad}
        #sens["con"] = {"x": [grad2]}

        return sens, fail


#--------------Make initial design variable dicts-----------------------#
DVdict = []
for dvname in ["area","aspect","camb0","cambf","ctwist", "dihedral","lesweep", "taper","tc0","tcf"]:
    tempDict = {"name" : dvname,
                "type" : "shape",
                "value" : 0.0,
                "capsGroup" : ""}
    DVdict.append(tempDict)

#setup thick DVs
capsGroup = ["rib","spar","OML"]
thickCt = 0
for dvname in ["thick1", "thick2", "thick3"]:
    tempDict = {"name" : dvname,
                "type" : "struct",
                "value" : 0.0,
                "capsGroup" : capsGroup[thickCt]}
    thickCt += 1
    DVdict.append(tempDict)

#----------------------Setup Optimization class---------------------------#
#initialize wing optimization class
wingOpt = WingOptimizer(DVdict)

# #test input::
# x = {}
# x["shape"] = [40.0, 6.0,  0.05, 0.05, 5.0,  5.0, 30.0,  0.5, 0.1, 0.1]
# x["struct"] = 0.01*np.ones(3)
# wingOpt.writeInput(x)

#setup pyOptSparse
sparseProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", wingOpt.objCon)

#shape DV settings
#["area","aspect","camb0","cambf","ctwist", "dihedral","lesweep", "taper","tc0","tcf"]
lbnds =   [20.0, 3.0, 0.01 , 0.01, 1.0,  1.0, 3.0,   0.3, 0.0, 0.0]
init = [40.0, 6.0,  0.05, 0.05, 5.0,  5.0, 30.0,  0.5, 0.1, 0.1]
ubnds =   [100.0,10.0, 0.3, 0.3, 10.0, 20.0, 50.0, 1.0, 0.3, 0.3]

#struct DV settings
#["rib","spar","OML"]
lBnds2 = 0.0001 * np.ones(3)
uBnds2 = 0.01*np.ones(3)
init2 = 0.001*np.ones(3)

#add design variable groups of each type
sparseProb.addVarGroup("shape", 10, lower=lbnds, upper=ubnds, value=init)
sparseProb.addVarGroup("struct", 3, lower=lBnds2, upper=uBnds2, value=init2)

#add functions, obj and constraint
sparseProb.addObj("obj")
#sparseProb.addConGroup("con", 1, lower=1, upper=1)

#setup SLSQP optimizer
optOptions = {"IPRINT": -1}
opt = SLSQP(options=optOptions)
sol = opt(sparseProb, sens=wingOpt.objGrad)

#print final solution to sol.out
wingOpt.write("Solution is...\n")
wingOpt.write("\txstar =  \n{}".format(sol.xStar))

wingOpt.close()