'''
script to run funtofem from design variable inputs, funtofem.input file and .csm files
returns function, gradient information in funtofem.out file
'''
from pyoptsparse import PSQP, Optimization
from caps2fun import Optimize as CapsOptimize
from caps2fun import Test
import numpy as np
import os

optimizationMode = "structural"

# #design parameters for airfoil
# despmtr diamAngle0 10
# despmtr chord0 1.0
# despmtr diamAnglef 10

# #design parameters of wing
# despmtr area 40.0
# despmtr taper 1.0

shapeActive = optimizationMode in ["TOGW"]
structActive = optimizationMode in ["structural","TOGW"]

#main script to run
    
DVnames = []

#--------------Make initial design variable dicts-----------------------#
DVdict = []
inits = [40.0, 1.0, 10.0, 10.0, 1.0]
ct = 0
DVind = 0
shapeind = 0
for dvname in ["area", "chord0","diamAngle0","diamAnglef","taper"]:
    tempDict = {"name" : dvname,
                "type" : "shape",
                "value" : inits[shapeind],
                "capsGroup" : "",
                "active" : shapeActive,
                "opt_ind" : DVind,
                "group_ind" : shapeind}
    DVdict.append(tempDict)
    if (shapeActive): DVind += 1
    shapeind += 1

##--------------setup thick DVs------------------##
#get csm file name from funtofem.cfg
cfgFile = os.path.join(os.getcwd(),"funtofem", "funtofem.cfg")
cfghdl = open(cfgFile, "r")
lines = cfghdl.readlines()
csmPrefix = ""
for line in lines:
    if ("csm" in line):
        chunks = line.split(" = ")
        csmPrefix = chunks[1].strip()
cfghdl.close()
csmFile = csmPrefix + ".csm"

#read csm file to count the current configuration
csmhdl = open(csmFile, "r")
lines = csmhdl.readlines()
csmhdl.close()
nribs = 0

for line in lines:
    chunks = line.split(" ")
    cfgpmtr = "cfgpmtr" in chunks[0]
    if (cfgpmtr):
        print(chunks)
        if ("nrib" in chunks[1]): nribs = int(chunks[2].strip())

nOML = nribs-1
nspars = 1

print(nribs,nspars,nOML)
def zeroString(nzeros):
    string = ""
    for i in range(nzeros):
        string += str(0)
    return string

groups = ["rib","spar","OML"]
numDVs = [nribs, nspars, nOML]
thickCt = 0
numThick = sum(numDVs)
numMaxDigits = len(str(numThick))

structind = 0

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

        #thickness = 0.001 * thickIndex #0.01
        initThickness = 0.01

        DVname = "thick" + zeroStr + str(thickIndex)
        DVnames.append(DVname)

        tempDict = {"name" : DVname,
                    "type" : "struct",
                    "value" : initThickness,
                    "capsGroup" : capsGroup,
                    "active" : structActive,
                    "opt_ind" : DVind,
                    "group_ind" : structind}
        if (structActive): DVind += 1
        thickCt += 1
        structind += 1
        DVdict.append(tempDict)

sorted = False

if (sorted):

    def capsCompare(elem):
        return elem["capsGroup"]


    structDVnames = []
    #make list of DV values and names
    for DV in DVdict:
        if (DV["type"] == "struct"):
            structDVnames.append(DV["name"])

    #sort the names and values based on ESP/CAPS sorting
    DVdict.sort(key=capsCompare)

    #get the sorted structDVs
    structDVs = []
    capsGroups = []
    for DV in DVdict:
        if (DV["type"] == "struct"):
            structDVs.append(DV["value"])
            capsGroups.append(DV["capsGroup"])

    #print(capsGroups)
    #print(structDVs)


#run a forward analysis to gauge the max stress constraint
mytest = Test(DVdict, functions=["ksfailure"])
functions = mytest.runForward()
maxStress = functions["ksfailure"]
maxStress = max(maxStress, 0.5)

#initialize wing optimization class
capsOpt = CapsOptimize(DVdict, optimizationMode)

#setup pyOptSparse
sparseProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", capsOpt.objCon)

#thickness in meters
names2 = DVnames
lbnds2 = 0.0001 * np.ones(thickCt)
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
sparseProb.addConGroup("con", 1,upper=maxStress)

#setup SLSQP optimizer
optOptions = {"IPRINT": -1}
opt = PSQP(options=optOptions)
sol = opt(sparseProb, sens=capsOpt.objGrad)


#after done with optimization
#print final solution to sol.out
capsOpt.cwrite("----------------------------")
capsOpt.cwrite("----------------------------\n")
capsOpt.cwrite("Opt Solution:  \n")
capsOpt.cwrite("\t{}".format(sol))

capsOpt.close()