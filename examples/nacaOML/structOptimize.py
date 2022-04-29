'''
script to run funtofem from design variable inputs, funtofem.input file and .csm files
returns function, gradient information in funtofem.out file
'''
from pyoptsparse import PSQP, Optimization
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

        #used DVs

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

        functions = ["ksfailure","mass"]
        for function in functions:
            line = "function,{}\n".format(function)
            inputHandle.write(line)

        #overwrite values in dict with the input x
        names = []
        values = []

        for DV in self.DVdict:
            #get attributes
            name = DV["name"]
            dvType = DV["type"]
            capsGroup = DV["capsGroup"]

            #get value
            if (DV["active"]):
                value = float(x[name])
                value = round(value, 5)

                #append values and names
                names.append(name)
                values.append(value)

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
        self.cwrite("----------------------------")
        self.cwrite("----------------------------\n")
        self.cwrite("Global Iteration #{}\n".format(self.iteration))
        self.cwrite("\tDV names {}\n".format(names))
        self.cwrite("\tDV values {}\n".format(values)) 
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
        #struct DV settings

    def readOutput(self, x):
        #read functions, gradients from funtofem call
        #make sure funtofem folder exists
        funtofemFolder = os.path.join(os.getcwd(), "funtofem")
        if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

        #make funtofem input file
        outputFile = os.path.join(funtofemFolder, "funtofem.out")
        outputHandle =  open(outputFile, "r")

        lines = outputHandle.readlines()
        

        self.functions = {}
        self.gradients = {}

        #read function values and gradients
        ifunc = -1
        firstLine = True
        for line in lines:
            if ("func" in line):
                #func, name, value
                ifunc += 1
                parts = line.split(",")
                name = parts[1]
                value = float(parts[2])

                self.functions[name] = value
                self.gradients[name] = np.zeros((nDV))
                iDV = 0

            elif ("grad" in line):
                #grad,name,deriv_i
                parts = line.split(",")
                name = parts[1]
                deriv = float(parts[2])
                self.gradients[name][iDV] = deriv
                iDV += 1

            elif (firstLine):
                firstLine = False
                #it's the first line
                parts = line.split(",")
                nfunc = int(parts[0])
                nDV = int(parts[1])

        #close the file
        outputHandle.close()

        #update status
        self.cwrite("\tRead output files for func, sens\n")

    def objCon(self, x):
        #py opt sparse function evaluator

        #write input, call funtofem, and get results
        self.writeInput(x)
        self.callFuntoFem()
        self.readOutput(x)

        #objective, constraint functions
        funcs = {}
        funcs["obj"] = self.functions["mass"]
        funcs["con"] = self.functions["ksfailure"]

        #write objective function
        self.cwrite("\tObjective func = {}\n".format(funcs["obj"]))

        fail = False

        return funcs, fail

    def objGrad(self, x, funcs):
        #pyoptsparse gradient function

        #read output just in case self.readOutput()
        objGrad = self.gradients["mass"]
        conGrad = self.gradients["ksfailure"]

        self.cwrite("\tObjective grad = {}\n".format(objGrad))

        #assume funtofem and esp/caps grid generation already done in objCon call
        fail = False

        sens = {}
        iDV = 0
        objSens = {}
        conSens = {}
        for DV in self.DVdict:
            name = DV["name"]
            if (name in self.usedDVs):
                objSens[name] = objGrad[iDV]
                conSens[name] = conGrad[iDV]
                iDV += 1

        sens["obj"] = objSens
        sens["con"] = conSens

        return sens, fail

    


#--------------Make initial design variable dicts-----------------------#
DVdict = []
inits = [40.0, 6.0,  0.05, 0.05, 5.0,  5.0, 0.0,  0.5, 0.1, 0.1]
ct = 0
for dvname in ["area","aspect","camb0","cambf","ctwist", "dihedral","lesweep", "taper","tc0","tcf"]:
    tempDict = {"name" : dvname,
                "type" : "shape",
                "value" : inits[ct],
                "capsGroup" : "",
                "active" : False}
    DVdict.append(tempDict)
    ct += 1

#setup thick DVs
groups = ["rib","spar","OML"]
numDVs = [8,2,14]
thickCt = 0
DVnames = []
for igroup in range(3):
    group = groups[igroup]
    numDV = numDVs[igroup]
    #print(group,numDV)

    for iDV in range(numDV):
        capsGroup = group + str(iDV+1)
        DVname = "thick" + str(thickCt + 1)
        DVnames.append(DVname)
        tempDict = {"name" : DVname,
                    "type" : "struct",
                    "value" : 0.001,
                    "capsGroup" : capsGroup,
                    "active" : True}
        thickCt += 1
        DVdict.append(tempDict)

#----------------------Setup Optimization class---------------------------#
#initialize wing optimization class
wingOpt = WingOptimizer(DVdict)
#print(DVdict)

#setup pyOptSparse
sparseProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", wingOpt.objCon)


#struct DV settings
#["rib","spar","OML"]
names2 = DVnames
lbnds2 = 0.0001 * np.ones(thickCt)
init2 = 0.001*np.ones(thickCt)
ubnds2 = 0.1*np.ones(thickCt)
scales2 = 100 * np.ones(thickCt)

for istruct in range(thickCt):
    name = names2[istruct]
    lbnd = lbnds2[istruct]
    init = init2[istruct]
    ubnd = ubnds2[istruct]
    scale = scales2[istruct]

    sparseProb.addVar(name, lower=lbnd, upper=ubnd, value = init, scale=scale)

#add functions, obj and constraint
sparseProb.addObj("obj")
#stress constraint upper bound 1/1.5
sparseProb.addConGroup("con", 1,upper=0.267)

#setup SLSQP optimizer
optOptions = {"IPRINT": -1}
opt = PSQP(options=optOptions)
sol = opt(sparseProb, sens=wingOpt.objGrad)


#after done with optimization
#print final solution to sol.out
wingOpt.cwrite("----------------------------")
wingOpt.cwrite("----------------------------\n")
wingOpt.cwrite("Opt Solution, xstar = \n")
wingOpt.cwrite("\t{}".format(D[:]))

wingOpt.close()