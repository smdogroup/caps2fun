import numpy as np
import matplotlib.pyplot as plt

aoa_file = "aoa_data.txt"
aoa_hdl = open(aoa_file, "r")
lines = aoa_hdl.readlines()
aoa_hdl.close()

functions = {}
aoa_vec = []
names = []
firstCase = False
storedAOA = False
storedData = False
for line in lines:
    chunks = line.split(",")
    chunk0 = chunks[0]
    if (len(chunks) > 1):
        chunk1 = chunks[1]
        chunk1_strip = chunk1.strip()
        if ("case" in line):
            caseind = int(chunk1_strip)
            if (firstCase): firstCase = False
            firstCase = caseind == 1
            if (storedAOA and not(storedData)): aoa_vec.pop()
            storedAOA = False
            storedData = False
        if ("func" in line):
            name = chunk1
            if (firstCase): 
                names.append(name)
                functions[name] = []
            value = chunks[2].strip()
            value = float(value)
            if (name == "temperature"):
                value = value + 300.0
                value /= 300.0
            elif (name == "mass"): 
                value /= 300.0
            functions[name].append(value)
            storedData = True
        if ("AOA" in line):
            aoa = float(chunk1_strip)
            aoa_vec.append(aoa)
            storedAOA = True

#reduce size of aoa_vec for empty ones
ndata = len(functions[name])
aoa_vec = aoa_vec[:ndata]
fig, ax = plt.subplots()
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')

print(aoa_vec)
colors = "bgr"
ct = 0
list1 = ["cl","cd","ksfailure"]
list2 = ["cl","cd"]
uselist1 = False
if (uselist1): 
    mylist = list1
else:
    mylist = list2
for name in mylist:
    color = colors[ct]; ct += 1
    if (name in "cl" + "cd"): 
        marker = "o"
        linestyle = "-"
    else:
        marker = "s"
        linestyle = ""
    plotArg = color + marker + linestyle
    #plot the curve
    ax.plot(aoa_vec, functions[name],plotArg, linewidth=3,label=name)
plt.xlabel("AOA (deg)")
plt.ylabel("analysis functions (ND)")
if (uselist1): plt.axis([-2.0, 10.0, 0, 0.6])
plt.legend()
plt.show()
plt.savefig("aoa_check.png")

