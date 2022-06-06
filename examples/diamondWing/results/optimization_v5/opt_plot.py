import numpy as np
import matplotlib.pyplot as plt

#Global Iteration #1
#	mass Obj = 3355.293036413158
#	ksfailure Con = 0.355864928437932

opt_hdl = open("opt_status.out","r")
lines = opt_hdl.readlines()
opt_hdl.close()

names = ["mass","ksfailure"]
iterations = []
functions = {}
functions["mass"] = []
functions["ksfailure"] = []
savedIteration = False
savedValues = False
for line in lines:
    if ("Global Iteration" in line):
        chunks = line.split("#")
        iteration = int(chunks[1].strip())
        iterations.append(iteration)
        savedIteration = True
    for name in names:
        if (name in line):
            chunks = line.split(" = ")
            value = float(chunks[1].strip())
            if (name == "mass"): value /= 1000.0
            functions[name].append(value)
            savedValues = True
iterations = iterations[:-1]

fix, ax = plt.subplots()
ax.plot
ax.axhline(y=0.267, color="g",linestyle="--", linewidth=2)
#ax.axvline(x=0, color='k')
colors = "bgcr"
ct = 0
for name in names:
    color = colors[ct]; ct += 1
    style = color + "o" + "-"
    ax.plot(iterations, functions[name], style, label=name, linewidth=3)
plt.xlabel("Optimization Iterations")
plt.ylabel("Normalized Functions")
plt.legend()
plt.show()
plt.savefig("sizingOpt_history.png")