from __future__ import print_function
import os

# ==============================================================================
# External Python modules
# ==============================================================================
from pprint import pprint

import numpy as np
from mpi4py import MPI

# ==============================================================================
# Extension modules
# ==============================================================================
from tacs import functions, constitutive, elements, TACS, pyTACS, problems

comm = MPI.COMM_WORLD

# Instantiate FEAAssembler
structOptions = {
    'printtimings':True,
    # Specify what type of elements we want in the f5
    #'writeSolution':True,
    #'outputElement': TACS.PLANE_STRESS_ELEMENT,
}

bdfFile = os.path.join(os.path.dirname(__file__), 'stiffPanel7.dat')
FEAAssembler = pyTACS(bdfFile, comm, options=structOptions)

def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    # Material properties
    rho = 2500.0        # density kg/m^3
    kappa = 230.0       # Thermal conductivity W/(m⋅K)
    specificHeat = 921.0 # Specific heat J/(kg⋅K)
    alpha = 8.6e-6       #CTE Titanium 6Al-4V
    E = 70e9            # Young's modulus (Pa)
    nu = 0.3            # Poisson's ratio
    ys = 464.0e6        # yield stress

    # Plate geometry
    tplate = 0.005    # 1 mm
    tMin = 0.0001    # 0.1 mm
    tMax = 0.05     # 5 cm

    # Setup property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, kappa=kappa, specific_heat=specificHeat, E=E, nu=nu, ys=ys)
    # Set one thickness dv for every component
    con = constitutive.IsoShellConstitutive(prop, t=tplate, tNum=dvNum, tlb=tMin, tub=tMax)
    #transform = elements.LinearElasticity3D(con)
    transform = None
    # For each element type in this component,
    # pass back the appropriate tacs element object
    elemList = []
    #model = elements.HeatConduction3D(con)
    for elemDescript in elemDescripts:
        if elemDescript in ['CQUAD4', 'CQUADR']:
            #basis = elements.LinearQuadBasis()
            elem = elements.Quad4NonlinearThermalShell(transform, con)
        elif elemDescript in ['CTRIA3', 'CTRIAR']:
            basis = elements.LinearTriangleBasis()
        else:
            print("Uh oh, '%s' not recognized" % (elemDescript))
        #elem = elements.Quad4ThermalShell(model, basis)
        elemList.append(elem)

    # Add scale for thickness dv
    scale = [100.0]
    return elemList, scale

# Set up constitutive objects and elements
FEAAssembler.initialize(elemCallBack)
transientOptions = {'printlevel':1}
# Setup problems
# Create a transient problem that will represent time varying convection
num_steps = 70
transientProb = FEAAssembler.createTransientProblem('Revised_NonlinThermalLoadAll_10imp_Transient', tInit=0.0, tFinal=35, numSteps=num_steps, options=transientOptions)

# Add problems to a list
allProblems = []

# Add functions to each problem
transientProb.addFunction('mass', functions.StructuralMass)
transientProb.addFunction('ks_temp', functions.KSTemperature,
                        ksWeight=100.0)


bdfInfo = FEAAssembler.getBDFInfo()
# cross-reference bdf object to use some of pynastrans advanced features
bdfInfo.cross_reference()
Pxy = []
eIDs = []
nIDs = []
fhz = 1.0
for nID in bdfInfo.nodes:
    nIDs.append(nID)

#nID_edge = [9,10,11,12,13,54,55,56,57,58,211,212,213,214,
#            238,239,240,241,526,527,528,529,553,554,555,
#            556,726,727,728,729,922,923,924,1094,1095,1096,1097]

nID_edge = [191,371,506,686,825,826,827,828,910,911,912,994,995,996,1078,1079,1080,1206,1210,1211,1212]


#Trying to add a slight initial geom imperfection

X = transientProb.getNodes()
Xpts = X
Lz = 1.0
Xpts[2::3] -= 0.01*np.sin(np.pi*Xpts[0::3]/Lz)
transientProb.setNodes(X)

#Adding dynamic loading
timeSteps = transientProb.getTimeSteps()
for step_i, time in enumerate(timeSteps):
    #F = np.array([0.0, 0.0, time*1000.0, 0.0, 0.0, 0.0]) #For Quad4Nonlinear/LinearShell
    F = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, time*2e-2]) #ForThermalBuckling deltaT=0.5s #steps=70 t_fin=35

    transientProb.addLoadToNodes(step_i, nIDs, F, nastranOrdering=True)
    #transientProb.addLoadToNodes(step_i, nID_edge, F, nastranOrdering=True)

allProblems.append(transientProb)

# Solve state for each problem, evaluate functions and sensitivities
funcs = {}
funcsSens = {}
for problem in allProblems:
    problem.solve()
    #problem.evalFunctions(funcs)
    #problem.evalFunctionsSens(funcsSens)
    problem.writeSolution()