import os
from caps2fun import Test, readnprocs

# #design parameters for airfoil
# despmtr diamAngle0 10
# despmtr chord0 1.0
# despmtr diamAnglef 10

# #design parameters of wing
# despmtr area 40.0
# despmtr taper 1.0

#main script to run
def makeDVdict():
    
    #settings
    shapeActive = False #def: False
    structActive = True #def: True

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
    DVnames = []
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
            thickness = 0.01

            DVname = "thick" + zeroStr + str(thickIndex)
            DVnames.append(DVname)

            tempDict = {"name" : DVname,
                        "type" : "struct",
                        "value" : thickness,
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

        print(capsGroups)
        print(structDVs)

    return DVdict

#-----------  Main Section, Run the Test -------------------#

#make a DV dict
DVdict = makeDVdict()

functions = ["ksfailure", "cl", "cd", "mass"]

#start a derivative test
mytest = Test(DVdict, functions=functions)
#mytest.runForward()
mytest.derivativeTest()