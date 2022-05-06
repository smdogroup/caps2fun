import os
import numpy as np

#read functions, gradients from funtofem call
#make sure funtofem folder exists
funtofemFolder = os.path.join(os.getcwd(), "funtofem")
if (not(os.path.exists(funtofemFolder))): os.mkdir(funtofemFolder)

#make funtofem input file
outputFile = os.path.join(funtofemFolder, "funtofem.out")
outputHandle =  open(outputFile, "r")

lines = outputHandle.readlines()


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
        gradients[functionName][iDV] = deriv
        iDV += 1

    elif (firstLine):
        firstLine = False
        #it's the first line
        parts = line.split(",")
        nfunc = int(parts[0])
        nDV = int(parts[1])

#close the file
outputHandle.close()

print("functions {}".format(functions))
print("gradients {}".format(gradients))