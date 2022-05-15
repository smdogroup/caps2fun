'''
script to run funtofem from design variable inputs, funtofem.input file and .csm files
returns function, gradient information in funtofem.out file
'''
import os, sys
import numpy as np

def writeInput(DVdict, functions, f2fanalysis, fun3danalysis, mode, eps=None, x_direction=None):
    # script to write the design variables to funtofem
    # sends them in the file funtofem.in
    # also tells which kind of analysis to perform, etc.

    #make sure funtofem folder exists
    funtofemFolder = os.path.join(os.getcwd(), "funtofem")
    if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

    #make funtofem input file
    inputFile = os.path.join(funtofemFolder, "funtofem.in")
    inputHandle =  open(inputFile, "w")

    #write the type of analysis to be done, "aerothermal","aerothermoelastic", "aeroelastic"
    #fun3d analysis : "inviscid", "laminar", "turbulent"
    line = "analysis,{},{}\n".format(f2fanalysis, fun3danalysis)
    inputHandle.write(line)

    #write the real or complex mode, aka adjoint or complex_step
    line = "mode,{}\n".format(mode)
    inputHandle.write(line)

    #write the analysis functions to be requested
    for function in functions:
        line = "function,{}\n".format(function)
        inputHandle.write(line)

    #overwrite values in dict with the input x
    names = []
    values = []

    nDV = 0
    for DV in DVdict:
        #get attributes
        name = DV["name"]
        dvType = DV["type"]
        capsGroup = DV["capsGroup"]

        #get value
        if (DV["active"]):
            nDV += 1

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


    #if mode is complex_step then also write the direction of complex step
    if (mode == "complex_step"):
        #if epsilon and direction are unspecified make default values
        if (eps is None): eps = 1.0e-30
        if (x_direction is None):
            x_direction = np.random.rand(nDV)
            x_direction = x_direction / np.linalg.norm(x_direction)

        #write the epsilon
        line = "epsilon,{}\n".format(eps)
        inputHandle.write(line)

        #write the direction vector
        inputHandle.write("x_dir,")
        for i in range(len(x_direction)):
            inputHandle.write("{}".format(x_direction[i]))
            if (i < len(x_direction)-1):
                inputHandle.write(",")

    #close the input file
    inputHandle.close()     

    
def makeDVdict():
    
    #settings
    shapeActive = False #def: False
    structActive = True #def: True

    #--------------Make initial design variable dicts-----------------------#
    DVdict = []
    inits = [40.0, 6.0,  0.05, 0.05, 5.0,  5.0, 0.0,  0.5, 0.1, 0.1]
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

    #print(DVdict)

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

def readOutput(DVdict,mode="adjoint"):
    
    #read functions, gradients from funtofem call
    #make sure funtofem folder exists
    funtofemFolder = os.path.join(os.getcwd(), "funtofem")
    if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

    #make funtofem input file
    outputFile = os.path.join(funtofemFolder, "funtofem.out")
    if (os.path.exists(outputFile)):
        #if the output file was written, the analysis ran successfully
        success = True

        #so read in the outputs
        outputHandle =  open(outputFile, "r")

        lines = outputHandle.readlines()
            
        if (mode == "adjoint"):

            functions = {}
            gradients = {}

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

                    functions[functionName] = value
                    gradients[functionName] = np.zeros((nDV))
                    iDV = 0

                elif ("grad" in line):
                    #grad,dvname,deriv_i
                    parts = line.split(",")
                    dvname = parts[1]
                    deriv = float(parts[2])

                    #find the DVind of that design variable (assuming out of order)
                    for DV in DVdict:
                        if (DV["name"] == dvname): ind = DV["opt_ind"]
                    #store the gradient at that index
                    gradients[functionName][ind] = deriv

                elif (firstLine):
                    firstLine = False
                    #it's the first line
                    parts = line.split(",")
                    nfunc = int(parts[0])
                    nDV = int(parts[1])

            #close the file
            outputHandle.close()

            #return the functions and gradients
            return functions, gradients

        elif (mode == "complex_step"):
            
            #read the funtofem.out with the following format
            #nfunc
            #func,fname,real,value,imag,value
            functions = {}

            for line in lines:
                chunks = line.split(",")
                if (len(chunks) == 1):
                    nfunc = chunks[0].strip()
                elif ("func" in line):
                    name = chunks[1]
                    real_part = float(chunks[3])
                    imag_part = float(chunks[5].strip())

                    functions[name] = real_part + imag_part * 1j

            #return the function values
            return functions


    else:
        sys.exit("Error: Funtofem Analysis failed, no funtofem.out file was created\n")