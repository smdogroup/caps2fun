'''
script to run funtofem from design variable inputs, funtofem.input file and .csm files
returns function, gradient information in funtofem.out file
'''
import os
import numpy as np

class Optimize():
    #class to run pyoptsparse optimize on the outside of funtofem

    def __init__(self, DVdict, optimizationMode):
        #set the DV dict here
        self.DVdict = DVdict

        #option: "structural", "full"
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
        #update design variable dict
        self.updateDVdict(x)

        if (self.optimizationMode == "structural"):
            funcs = ["ksfailure", "mass"]
        elif (self.optimizationMode == "full"):
            funcs = ["ksfailure","cl","cd","mass"]

        writeInput(self.DVdict, funcs)

        #update status to the inputs that were run
        self.cwrite("----------------------------")
        self.cwrite("----------------------------\n")
        self.cwrite("Global Iteration #{}\n".format(self.iteration))
        self.cwrite("\tDV names {}\n".format(names))
        self.cwrite("\tDV values {}\n".format(values))              

    def updateDVdict(self,x):
        tempDict = []
        for DV in self.DVdict:

            #overwrite the value
            if (DV["active"]):
                name = DV["name"]
                DV["value"] = float(x[name])

            tempDict.append(DV)
        
        #overwrite DVdict
        self.DVdict = tempDict

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

        #self.deleteF2Ffiles()

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