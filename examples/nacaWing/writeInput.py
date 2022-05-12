'''
script to run funtofem from design variable inputs, funtofem.input file and .csm files
returns function, gradient information in funtofem.out file
'''
import os
import numpy as np

def writeInput(DVdict, functions):
    # script to write the design variables to funtofem
    # sends them in the file funtofem.in
    # also tells which kind of analysis to perform, etc.

    #make sure funtofem folder exists
    funtofemFolder = os.path.join(os.getcwd(), "funtofem")
    if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

    #make funtofem input file
    inputFile = os.path.join(funtofemFolder, "funtofem.in")
    inputHandle =  open(inputFile, "w")

    for function in functions:
        line = "function,{}\n".format(function)
        inputHandle.write(line)

    #overwrite values in dict with the input x
    names = []
    values = []

    for DV in DVdict:
        #get attributes
        name = DV["name"]
        dvType = DV["type"]
        capsGroup = DV["capsGroup"]

        #get value
        if (DV["active"]):
            value = DV["value"]
            #value = round(value, 5)

            #append values and names
            names.append(name)
            values.append(value)

            #make the line strings
            if (dvType == "shape"):
                line = "{},{},{:5f}\n".format(name, dvType, value)
            elif (dvType == "struct"):
                line = "{},{},{},{:5f}\n".format(name, dvType, capsGroup, value)   

            #write the line to the input file
            inputHandle.write(line)


    #close the input file
    inputHandle.close()     

    


#--------------- Settings --------------------------------------------#
shapeActive = False #def: False
structActive = True #def: True


#--------------Make initial design variable dicts-----------------------#
DVdict = []
inits = [40.0, 6.0,  0.05, 0.05, 5.0,  5.0, 0.0,  0.5, 0.1, 0.1]
ct = 0
DVind = 0
for dvname in ["area","aspect","camb0","cambf","ctwist", "dihedral","lesweep", "taper","tc0","tcf"]:
    tempDict = {"name" : dvname,
                "type" : "shape",
                "value" : inits[ct],
                "capsGroup" : "",
                "active" : shapeActive,
                "ind" : DVind}
    DVdict.append(tempDict)
    DVind += 1
    ct += 1

#setup thick DVs
nribs = 16
nspars = 2
nOML = nribs-1

def zeroString(nzeros):
    string = ""
    for i in range(nzeros):
        string += str(0)
    return string

groups = ["rib","spar","OML"]
numDVs = [nribs, nspars, nOML]
thickCt = 0
DVnames = []
numThick = sum(numDVs)
numMaxDigits = len(str(numThick))

DVind = 0

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

        thickness = 0.001 * thickIndex #0.01

        DVname = "thick" + zeroStr + str(thickIndex)
        DVnames.append(DVname)

        tempDict = {"name" : DVname,
                    "type" : "struct",
                    "value" : thickness,
                    "capsGroup" : capsGroup,
                    "active" : structActive,
                    "ind" : DVind}
        DVind += 1
        thickCt += 1
        DVdict.append(tempDict)

functions = ["ksfailure","cl","cd","mass"]

print(DVdict)

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

    print(capsGroups)
    print(structDVs)

writeInput(DVdict, functions)