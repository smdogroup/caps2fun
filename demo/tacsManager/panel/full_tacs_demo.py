"""
full verbatim demo of what pytacs is doing here
"""

from tacs import functions, constitutive, elements, pyTACS
from pyNastran.bdf.bdf import read_bdf
import caps2tacs
from mpi4py import MPI

# all the capsWrapper package modules/classes are also available 
# in caps2tacs since caps2tacs package imports the whole capsWrapper package\

panel_problem = caps2tacs.CapsStruct.default(csmFile="panel.csm")
tacs_aim = panel_problem.tacsAim
egads_aim = panel_problem.egadsAim

madeupium = caps2tacs.Isotropic.madeupium()
tacs_aim.add_material(material=madeupium)


constraint = caps2tacs.ZeroConstraint(name="fixEdge", caps_constraint="edge")
tacs_aim.add_constraint(constraint=constraint)

load = caps2tacs.GridForce(name="load1", caps_load="plate", direction=[0,0,-1.], magnitude=1.0E3)    
tacs_aim.add_load(load=load)

thick_DV = caps2tacs.ThicknessVariable(name="thick", caps_group="plate", value=0.01, material=madeupium)
tacs_aim.add_variable(variable=thick_DV)

despmtrs = ["plateLength", "plateWidth"]
for despmtr in despmtrs:
    shape_var = caps2tacs.ShapeVariable(name=despmtr)
    tacs_aim.add_variable(variable=shape_var)

tacs_aim.setup_aim()
egads_aim.set_mesh()

# start a caps tacs main problem
caps_tacs = caps2tacs.CapsTacs(tacs_aim=tacs_aim, egads_aim=egads_aim)

"""
next section is full description of what pytacs does and we will try to simplify it
"""

#Step 1) initialize pytacs with that data file
#FEASolver = pyTACS(caps_tacs.dat_file)
dat_file = caps_tacs.dat_file
comm = MPI.COMM_WORLD

# Step 2) open up pyMeshLoader
# self.meshLoader = pyMeshLoader(self.comm, self.dtype, debugFlag)
bdfInfo = read_bdf(dat_file, validate=False, xref=False, debug=False)

# self.meshLoader.scanBdfFile(fileName)
elemConnectivity = [None] * bdfInfo.nelements
elemConnectivityPointer = [None] * (bdfInfo.nelements + 1)
elemConnectivityPointer[0] = 0
elementObjectCounter = 0
# List specifying which tacs element object each element in bdf should point to
elemObjectNumByElem = [None] * (bdfInfo.nelements)

bdfXpts = bdfInfo.get_xyz_in_coord()

# Loop through every element and record information needed for tacs
for tacsElementID, nastranElementID in enumerate(bdfInfo.element_ids):
    element = bdfInfo.elements[nastranElementID]
    elementType = element.type.upper()
    propertyID = element.pid
    componentID = self.idMap(propertyID, self.nastranToTACSCompIDDict)

    # This element type has not been added to the list for the component group yet, so we append it
    if elementType not in self.elemDescripts[componentID]:
        self.elemDescripts[componentID].append(elementType)
        self.elemObjectNumByComp[componentID].append(elementObjectCounter)
        elementObjectCounter += 1

    # Find the index number corresponding to the element object number for this component
    componentTypeIndex = self.elemDescripts[componentID].index(elementType)
    self.elemObjectNumByElem[tacsElementID] = self.elemObjectNumByComp[componentID][componentTypeIndex]

    # We've identified a ICEM property label
    if 'Shell element data for family' in element.comment:
        componentName = element.comment.split()[-1]
        self.compDescripts[componentID] = componentName

    conn = element.nodes.copy()

    # TACS has a different node ordering than Nastran for certain elements,
    # we now perform the reordering (if necessary)
    if elementType in ['CQUAD4', 'CQUADR']:
        conn = [conn[0], conn[1], conn[3], conn[2]]
    elif elementType in ['CQUAD9', 'CQUAD']:
        conn = [conn[0], conn[4], conn[1], conn[7], conn[8], conn[5], conn[3], conn[6], conn[2]]
    elif elementType in ['CHEXA8', 'CHEXA']:
        conn = [conn[0], conn[1], conn[3], conn[2], conn[4], conn[5], conn[7], conn[6]]

    # Map node ids in connectivity from Nastan numbering to TACS numbering
    self.elemConnectivity[tacsElementID] = self.idMap(conn, self.nastranToTACSNodeIDDict)
    self.elemConnectivityPointer[tacsElementID + 1] = self.elemConnectivityPointer[
                                                            tacsElementID] + len(element.nodes)

# Allocate list for user-specified tacs element objects
self.elemObjects = [None] * elementObjectCounter


# # Save pynastran bdf object
# self.bdfInfo = self.meshLoader.getBDFInfo()

# Set up TACS Assembler
#FEASolver.initialize()

#choose the functions to evaluate
#evalFuncs = ['wing_mass', 'ks_vmfailure']

#read the bdf & dat file into pytacs FEAsolver
#SPs represents "StructuralProblems"
#SPs = FEASolver.createTACSProbsFromBDF()

# Read in forces from BDF and create tacs struct problems
# for caseID in SPs:
#     SPs[caseID].addFunction('wing_mass', functions.StructuralMass)
#     SPs[caseID].addFunction('ks_vmfailure', functions.KSFailure, safetyFactor=1.5, ksWeight=1000.0)

# Solve each structural problem and write solutions
# func = {}; sens = {}
# for caseID in SPs:
#     SPs[caseID].solve()
#     SPs[caseID].evalFunctions(func,evalFuncs=evalFuncs)
#     SPs[caseID].evalFunctionsSens(sens,evalFuncs=evalFuncs)
#     SPs[caseID].writeSolution(outputDir=os.path.dirname(__file__))