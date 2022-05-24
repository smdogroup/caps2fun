'''
script to run funtofem from design variable inputs, funtofem.input file and .csm files
returns function, gradient information in funtofem.out file
'''
from pyoptsparse import PSQP, Optimization
from caps2fun import Optimize as CapsOptimize
import numpy as np
import os

optimizationMode = "structural"

## ------------------ Make DV Dict --------------------------------- ##

shapeActive = optimizationMode == "full"
structActive = True

initThickness = 0.10; #10 mm
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

#initialize wing optimization class
capsOpt = CapsOptimize(DVdict, optimizationMode)

#setup pyOptSparse
sparseProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", capsOpt.objCon)

#thickness in meters
names2 = DVnames
lbnds2 = 0.001 * np.ones(thickCt)
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
sparseProb.addConGroup("con", 1,upper=0.10)

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