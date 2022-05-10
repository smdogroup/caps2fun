import pyCAPS
import os
import argparse

csmFile = "naca_OML_fluid.csm"
caps = pyCAPS.Problem(problemName="fluid",
                    capsFile=csmFile,
                    outLevel=1)

#load meshing aims
caps.analysis.create(aim = "egadsTessAIM",name="egadsTess")
meshAIM = caps.analysis.create(aim = "tetgenAIM", name = "tetgen")

#set new egads body tesselation parameters
caps.analysis["egadsTess"].input.Tess_Params = [15.0, 0.01, 20.0]
meshAIM.input.Preserve_Surf_Mesh = True
meshAIM.input["Surface_Mesh"].link(caps.analysis["egadsTess"].output["Surface_Mesh"])
meshAIM.input.Mesh_Format = "AFLR3"

fun3dAIM = caps.analysis.create(aim = "fun3dAIM",
                                name = "fun3d")
# Set project name
fun3dAIM.input.Proj_Name = "fun3dTetgenTest"

# Link the mesh
fun3dAIM.input["Mesh"].link(meshAIM.output["Volume_Mesh"])

fun3dAIM.input.Mesh_ASCII_Flag = False

# Set AoA number
caps.analysis["fun3d"].input.Alpha = 0.0

# Set Mach number
caps.analysis["fun3d"].input.Mach = 0.7

# Set equation type
fun3dAIM.input.Equation_Type = "incompressible"

# Set Viscous term
caps.analysis["fun3d"].input.Viscous = "inviscid"

# Set number of iterations
caps.analysis["fun3d"].input.Num_Iter = 10

# Set CFL number schedule
caps.analysis["fun3d"].input.CFL_Schedule = [0.5, 3.0]

# Set read restart option
fun3dAIM.input.Restart_Read = "off"

# Set CFL number iteration schedule
caps.analysis["fun3d"].input.CFL_Schedule_Iter = [1, 40]

# Set overwrite fun3d.nml if not linking to Python library
caps.analysis["fun3d"].input.Overwrite_NML = True


inviscidBC1 = {"bcType" : "Inviscid", "wallTemperature" : 1}
fun3dAIM.input.Boundary_Condition = {"OML"   : inviscidBC1,
                                   "Symmetry" : "SymmetryY",
                                     "Farfield": "Farfield"}

fun3dAIM.preAnalysis()

print("\n\nRunning FUN3D.....")
currentDirectory = os.getcwd()

os.chdir(fun3dAIM.analysisDir)
os.system("nodet_mpi --animation_freq -1 --write_aero_loads_to_file> Info.out")
os.chdir(currentDirectory)

fun3dAIM.postAnalysis()

print ("Total Force - Pressure + Viscous")
# Get Lift and Drag coefficients
print ("Cl = " , fun3dAIM.output.CLtot,
       "Cd = " , fun3dAIM.output.CDtot)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx = " , fun3dAIM.output.CMXtot,
       "Cmy = " , fun3dAIM.output.CMYtot,
       "Cmz = " , fun3dAIM.output.CMZtot)

# Get Cx, Cy, Cz coefficients
print ("Cx = " , fun3dAIM.output.CXtot,
       "Cy = " , fun3dAIM.output.CYtot,
       "Cz = " , fun3dAIM.output.CZtot)

print ("Pressure Contribution")
# Get Lift and Drag coefficients
print ("Cl_p = " , fun3dAIM.output.CLtot_p,
       "Cd_p = " , fun3dAIM.output.CDtot_p)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx_p = " , fun3dAIM.output.CMXtot_p,
       "Cmy_p = " , fun3dAIM.output.CMYtot_p,
       "Cmz_p = " , fun3dAIM.output.CMZtot_p)

# Get Cx, Cy, and Cz, coefficients
print ("Cx_p = " , fun3dAIM.output.CXtot_p,
       "Cy_p = " , fun3dAIM.output.CYtot_p,
       "Cz_p = " , fun3dAIM.output.CZtot_p)

print ("Viscous Contribution")
# Get Lift and Drag coefficients
print ("Cl_v = " , fun3dAIM.output.CLtot_v,
       "Cd_v = " , fun3dAIM.output.CDtot_v)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx_v = " , fun3dAIM.output.CMXtot_v,
       "Cmy_v = " , fun3dAIM.output.CMYtot_v,
       "Cmz_v = " , fun3dAIM.output.CMZtot_v)

# Get Cx, Cy, and Cz, coefficients
print ("Cx_v = " , fun3dAIM.output.CXtot_v,
       "Cy_v = " , fun3dAIM.output.CYtot_v,
       "Cz_v = " , fun3dAIM.output.CZtot_v)
