import matplotlib.pyplot as plt
import numpy as np

nelems = np.array([4669, 16874, 75667])
#log_nelems = np.log(nelems)

filename = "conv_tess1.txt" #"convergence_stats.txt"
conv_hdl = open(filename, "r")
lines = conv_hdl.readlines()
conv_hdl.close()

cases = []
mean_values = {}
stddev_values = {}
mesh_values = {}
inCases = False
for line in lines:
    chunks = line.split(",")
    chunk0 = chunks[0]
    print(line)
    if (len(chunks)>1):
        chunk1 = chunks[1].strip()
        if ("case" in chunk0):
            inCases = True
            caseind = int(chunk1)
            caseDict = {}
            caseDict["ind"] = caseind
        elif ("func" in chunk0):
            name = chunks[1]
            if (len(chunks) > 4): 
                mean = float(chunks[3])
                stddev = float(chunks[5])
            else:
                mean = float(chunks[2])
                stddev = 0.0
            if (name == "temperature"):
                mean /= 5.0 * 300
                stddev /= 5.0 * 300
            if (caseind == 1):
                mean_values[name] = []
                stddev_values[name] = []
            mean_values[name].append(mean)
            stddev_values[name].append(stddev)
            caseDict["mean"] = mean
            caseDict["stddev"] = stddev
        elif ("nelems" in chunk0):
            nelems = int(chunk1)
            if (caseind == 1):
                mesh_values["nelems"] = []
            mesh_values["nelems"].append(nelems)
            caseDict["nelems"] = nelems
        elif ("tess1" in chunk0):
            tess1 = float(chunk1)
            if (caseind == 1):
                mesh_values["tess1"] = []
            mesh_values["tess1"].append(tess1)
            caseDict["tess1"] = nelems
    else:
        if (inCases):
            inCases = False
            cases.append(caseDict)

print(cases)       

colors = "kbgr"
colorind = 0
mesh_string = "nelems" #nelems, tess1
for name in ["ksfailure","cl","cd","temperature"]:
    color = colors[colorind]; colorind += 1
    if (name in ["ksfailure","temperature"]):
        linestyle = "-"
        marker = "o"
    else:
        linestyle = "--"
        marker = "s"
    if (filename == "convergence_stats.txt"): linestyle = ""
    stylestr = color + marker + linestyle
    plt.plot(mesh_values[mesh_string], mean_values[name], stylestr, label = name, linewidth=3)
    #print(mean_values[name],stddev_values[name])
    plt.errorbar(mesh_values[mesh_string], mean_values[name], yerr=stddev_values[name], fmt=marker,ecolor = 'red',color=color)

plt.xlabel(mesh_string + " struct mesh")
plt.ylabel("normalized function values")
if (mesh_string == "nelems"): plt.axis([1e3, 1e5, -2, 6])
plt.xscale('log')
plt.legend()
plt.show()
plt.savefig("struct_conv1.png")