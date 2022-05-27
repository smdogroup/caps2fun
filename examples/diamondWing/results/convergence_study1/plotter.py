import matplotlib.pyplot as plt
import numpy as np

nelems = np.array([4669, 16874, 75667])
#log_nelems = np.log(nelems)

conv_hdl = open("convergence_stats.txt", "r")
lines = conv_hdl.readlines()
conv_hdl.close()

cases = []
mean_values = {}
stddev_values = {}
nelem_values = []
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
            mean = float(chunks[3])
            stddev = float(chunks[5])
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
            nelem_values.append(nelems)
            caseDict["nelems"] = nelems
    else:
        if (inCases):
            inCases = False
            cases.append(caseDict)

print(cases)       

colors = "kbcg"
colorind = 0
for name in ["ksfailure","cl","cd","temperature"]:
    color = colors[colorind]; colorind += 1
    if (name in ["ksfailure","temperature"]):
        linestyle = "-"
    else:
        linestyle = "--"
    stylestr = color + "o" + linestyle
    plt.plot(nelem_values, mean_values[name], stylestr, label = name, linewidth=3)
    #print(mean_values[name],stddev_values[name])
    plt.errorbar(nelem_values, mean_values[name], yerr=stddev_values[name], fmt='o',ecolor = 'red',color=color)

plt.xlabel("nelems struct mesh")
plt.ylabel("normalized function values")
plt.axis([1e3, 1e5, -2, 6])
plt.xscale('log')
plt.legend()
plt.show()
plt.savefig("struct_conv1.png")