import os
from caps2fun import Test

#main script to run
def makeDVdict():
    
    #settings
    shapeActive = False #def: False
    structActive = True #def: True

    #--------------Make initial design variable dicts-----------------------#
    DVdict = []
    inits = [120.0, 6.0,  0.03, 0.03, 5.0,  5.0, 0.0,  0.5, 0.05, 0.05]
    ct = 0
    DVind = 0
    shapeind = 0
    for dvname in ["area","aspect","camb0","cambf","ctwist", "dihedral","lesweep", "taper","tc0","tcf"]:
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
    cfgFile = os.path.join(os.getcwd(),"caps2fun.cfg")
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
    nspars = 0
    stringerOn = 0

    for line in lines:
        chunks = line.split(" ")
        cfgpmtr = "cfgpmtr" in chunks[0]
        if (cfgpmtr):
            print(chunks)
            if ("nrib" in chunks[1]): nribs = int(chunks[2].strip())
            if ("nspar" in chunks[1]): nspars = int(chunks[2].strip())
            if ("stringerOn" in chunks[1]): stringerOn = int(chunks[2].strip())

    nOML = nribs-1

    print(nribs,nspars,stringerOn,nOML)

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
            thickness = 0.010

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

    
    #add stringer caps group DV if turned on
    if (stringerOn == 1):

        thickIndex = thickCt + 1
        numDigits = len(str(thickIndex))
        numZeros = numMaxDigits - numDigits
        zeroStr = zeroString(numZeros)

        #thickness = 0.001 * thickIndex #0.01
        thickness = 0.010

        DVname = "thick" + zeroStr + str(thickIndex)
        DVnames.append(DVname)

        #tempDict for stringers
        tempDict = {"name" : DVname,
                    "type" : "struct",
                    "value" : thickness,
                    "capsGroup" : "stringer",
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

functions = ["ksfailure", "temperature","cl", "cd", "mass"]

# #start a derivative test
mytest = Test(DVdict, functionNames=functions)
mytest.runForward()

