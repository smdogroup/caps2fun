from pyoptsparse import SLSQP, Optimization
import numpy as np

def objCon(x):
    fail = False
    x1 = x["struct"]
    value = x1[0]**2 + x1[1]**2 + x1[2]**2
    x2 = x["shape"]
    value += x2[0]**2 + x2[1]**2 + x2[2]**2
    funcs = {}
    funcs["obj"] = value
    return funcs, fail
def objGrad(x, funcs):
    x1 = x["struct"]
    x2 = x["shape"]
    sens = {}
    grad = 2*x1
    grad2 = 2*x2
    #grad2 = 0*x2
    sens = {"obj" : {"struct" : grad, "shape" : grad2}}
    print(grad)
    return sens
#setup pyOptSparse
sparseProb = Optimization("Stiffened Panel Aerothermoelastic Optimization", objCon)

#  ["area","aspect","camb0","cambf","ctwist", "dihedral","lesweep", "taper","tc0","tcf"]
lbnds =   [20.0, 3.0, 0.0 , 0.0, 1.0,  1.0, 3.0,   0.3, 0.0, 0.0]
init = [40.0, 6.0,  0.0, 0.0, 5.0,  5.0, 30.0,  0.5, 0.1, 0.1]
ubnds =   [100.0,10.0, 0.3, 0.3, 10.0, 20.0, 50.0, 1.0, 0.3, 0.3]

#["rib","spar","OML"]
lBnds2 = 0.0001 * np.ones(3)
uBnds2 = 0.01*np.ones(3)
init2 = 0.001*np.ones(3)

#sparseProb.addVarGroup("shape", 10, "c", lower=lbnds, upper=ubnds, value=init)
sparseProb.addVarGroup("struct", 3,lower=lBnds2, upper=uBnds2, value=init2)
sparseProb.addVarGroup("shape", 3, lower=lBnds2, upper=uBnds2, value=init2)

#sparseProb.addConGroup("con", 1, lower=1, upper=1)
sparseProb.addObj("obj")

# if comm.rank == 0:
#     print(sparseProb)

#optOptions = {"IPRINT": -1}
opt = SLSQP(options={})
sol = opt(sparseProb, sens=objGrad)

print(sol)
print('\nsol.xStar:  ', sol.xStar)