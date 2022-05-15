'''
script to run funtofem from design variable inputs, funtofem.input file and .csm files
returns function, gradient information in funtofem.out file
'''
from pyoptsparse import PSQP, Optimization
import os
import numpy as np

class WingOptimizer():
    #class to run pyoptsparse optimize on the outside of funtofem

    def __init__(self, DVdict, optimizationMode):
        #set the DV dict here
        self.DVdict = DVdict
        self.optimizationMode = optimizationMode

        #other information such as analysis type, mesher type is  setup in funtofem.py
        self.n_procs = 192

        #make status file
        statusFile = os.path.join(os.getcwd(), "opt_status.out")
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
        self.inputFile = os.path.join(funtofemFolder, "funtofem.in")
        inputHandle =  open(self.inputFile, "w")

        #write the type of analysis to be done, in this case adjoint
        line = "mode,{}\n".format("adjoint")

        #write the analysis functions to be performed
        if (self.optimizationMode == "structural"):
            functions = ["ksfailure","mass"]
        elif (self.optimizationMode == "full"):
            functions = ["cl","cd","ksfailure","mass"]

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
                #value = round(value, 5)

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

    def readOutput(self):
        #read functions, gradients from funtofem call
        #make sure funtofem folder exists
        funtofemFolder = os.path.join(os.getcwd(), "funtofem")
        if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

        #make funtofem input file
        self.outputFile = os.path.join(funtofemFolder, "funtofem.out")
        if (os.path.exists(self.outputFile)):
            #if the output file was written, the analysis ran successfully
            self.success = True

            #so read in the outputs
            outputHandle =  open(self.outputFile, "r")

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
                    functionName = parts[1]
                    value = float(parts[2])

                    self.functions[functionName] = value
                    self.gradients[functionName] = np.zeros((nDV))
                    iDV = 0

                elif ("grad" in line):
                    #grad,dvname,deriv_i
                    parts = line.split(",")
                    dvname = parts[1]
                    deriv = float(parts[2])

                    #find the DVind of that design variable (assuming out of order)
                    for DV in self.DVdict:
                        if (DV["name"] == dvname):
                            ind = DV["opt_ind"]
                    self.gradients[functionName][ind] = deriv
                    iDV += 1

                elif (firstLine):
                    firstLine = False
                    #it's the first line
                    parts = line.split(",")
                    nfunc = int(parts[0])
                    nDV = int(parts[1])

            #close the file
            outputHandle.close()

            #delete output file
            os.remove(self.outputFile)

            #update status
            self.cwrite("\tRead output files for func, sens\n")

            #since didn't fail update iteration count
            self.iteration += 1
            self.fail = 0

        else:
            #otherwise if the output file was not written, it failed
            #possibly due to negative volume or something
            #can try and change parameters here or pass information to change settings and rerun
            self.success = False
            self.cwrite("\tAnalysis Failed\n")
            self.fail = 1

        return self.fail

    def deleteF2Ffiles(self):
        os.remove(self.inputFile)
        os.remove(self.outputFile)

    def objCon(self, x):
        #py opt sparse function evaluator

        #write input, call funtofem, and get results
        self.writeInput(x)
        self.callFuntoFem()
        self.readOutput()

        self.deleteF2Ffiles()

        #objective, constraint functions
        funcs = {}

        if (self.optimizationMode == "structural"):
            funcs["obj"] = self.functions["mass"]
            funcs["con"] = self.functions["ksfailure"]
            objname = "mass"
            conname = "ksfailure"
        elif (self.optimizationMode == "full"):
            #TBD some kind of min fuel burn thing
            funcs["obj"] = 1
            funcs["con"] = 1
            objname = ""
            conname = ""

        #write objective function
        self.cwrite("\t{} Obj = {}\n".format(objname,funcs["obj"]))
        self.cwrite("\t{} Con = {}\n".format(conname,funcs["con"]))

        return funcs, self.fail

    def objGrad(self, x, funcs):
        #pyoptsparse gradient function

        if (self.optimizationMode == "structural"):
            objGrad = self.gradients["mass"]
            conGrad = self.gradients["ksfailure"]
            objname = "mass"
            conname = "ksfailure"
        elif (self.optimizationMode == "full"):
            #TBD some kind of min fuel burn thing
            #objGrad = self.gradients["mass"]
            #conGrad = self.gradients["ksfailure"]
            objname = ""
            conname = ""
        

        self.cwrite("\t {} Obj grad = {}\n".format(objname,objGrad))
        self.cwrite("\t{} Con grad = {}\n".format(conname,conGrad))

        sens = {}
        iDV = 0
        objSens = {}
        conSens = {}
        for DV in self.DVdict:
            name = DV["name"]
            if (DV["active"]):
                objSens[name] = objGrad[iDV]
                conSens[name] = conGrad[iDV]
                iDV += 1

        sens["obj"] = objSens
        sens["con"] = conSens

        return sens, self.fail



#--------------------  Settings ----------------------------------------#
#structural optimization "structural" or full optimization "full"
optimizationMode = "structural" #"full"

#--------------Make initial design variable dicts-----------------------#
shapeActive = False
initThickness = 0.01; #10 mm

if (optimizationMode == "full"): shapeActive = True
DVdict = []
inits = [40.0, 6.0,  0.05, 0.05, 5.0,  5.0, 0.0,  0.5, 0.1, 0.1]
ct = 0
DVind = 0
ishape = 0
istruct = 0
for dvname in ["area","aspect","camb0","cambf","ctwist", "dihedral","lesweep", "taper","tc0","tcf"]:
    tempDict = {"name" : dvname,
                "type" : "shape",
                "value" : inits[ct],
                "capsGroup" : "",
                "active" : shapeActive,
                "opt_ind" : DVind,
                "group_ind" : ishape}
    if (shapeActive): DVind += 1
    ishape += 1
    DVdict.append(tempDict)

    ct += 1

#setup thick DVs
structActive = False
if (optimizationMode in ["structural","full"]): structActive = True

def zeroString(nzeros):
    string = ""
    for i in range(nzeros):
        string += str(0)
    return string
 
#read the number of ribs and spars from the file
f = open("nacaWing.csm","r")
lines = f.readlines()
nribs = 0
nspars = 0
for line in lines:
    isConfig = "cfgpmtr" in line
    chunks = line.split(" ")
    if (isConfig and "nrib" in line):
        nribs = int(chunks[-1].strip())
    if (isConfig and "nspar" in line):
        nspars = int(chunks[-1].strip())
f.close()
nOML = nribs-1

groups = ["rib","spar","OML"]
numDVs = [nribs, nspars, nOML]
thickCt = 0
DVnames = []
numThick = sum(numDVs)
numMaxDigits = len(str(numThick))

for igroup in range(3):
    group = groups[igroup]
    numDV = numDVs[igroup]
    #print(group,numDV)

    for iDV in range(numDV):
        capsGroup = group + str(iDV+1)
        thickIndex = thickCt + 1
        numDigits = len(str(thickIndex))
        numZeros = numMaxDigits - numDigits
        zeroStr = zeroString(numZeros)
        DVname = "thick" + zeroStr + str(thickIndex)
        DVnames.append(DVname)
        tempDict = {"name" : DVname,
                    "type" : "struct",
                    "value" : initThickness,
                    "capsGroup" : capsGroup,
                    "active" : structActive,
                    "opt_ind" : DVind,
                    "group_ind" : istruct}
        if (structActive): DVind += 1
        thickCt += 1
        istruct += 1
        DVdict.append(tempDict)

#end of design variable initialization section

#---------------------- Run Optimization ---------------------------#

#initialize wing optimization class
wingOpt = WingOptimizer(DVdict, optimizationMode)

#setup pyOptSparse
sparseProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", wingOpt.objCon)


#struct DV settings
#["rib","spar","OML"]

#thickness in meters
names2 = DVnames
lbnds2 = 0.0005 * np.ones(thickCt)
init2 = initThickness*np.ones(thickCt)
ubnds2 = 1.0*np.ones(thickCt)
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
#stress constraint upper bound 1/1.5/2.5 = 0.267
sparseProb.addConGroup("con", 1,upper=0.267)

#setup SLSQP optimizer
optOptions = {"IPRINT": -1}
opt = PSQP(options=optOptions)
sol = opt(sparseProb, sens=wingOpt.objGrad)


#after done with optimization
#print final solution to sol.out
wingOpt.cwrite("----------------------------")
wingOpt.cwrite("----------------------------\n")
wingOpt.cwrite("Opt Solution:  \n")
wingOpt.cwrite("\t{}".format(sol))

wingOpt.close()