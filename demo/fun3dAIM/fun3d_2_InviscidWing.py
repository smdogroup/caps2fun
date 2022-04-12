#------------------------------------------------------------------------------#

# Import pyCAPS module
import pyCAPS

# Import os module
import os

# f90nml is used to write fun3d inputs not available in the aim
import f90nml

#------------------------------------------------------------------------------#

# Load CSM file
filename = os.path.join("naca_small.csm")
myProblem = pyCAPS.Problem(problemName = "workDir_01_Defaults",
                           capsFile = filename,
                           outLevel = 1)

# Alias the geometry
transport = myProblem.geometry

# Change to Inviscid CFD view
transport.cfgpmtr.cfdOn = 1

#------------------------------------------------------------------------------#

# Create aflr4 AIM
aflr4 = myProblem.analysis.create(aim  = "aflr4AIM",
                                  name = "aflr4")

#aflr4.geometry.view()

# Farfield growth factor
aflr4.input.ff_cdfr = 1.4

# Scaling factor to compute AFLR4 'ref_len' parameter
aflr4.input.Mesh_Length_Factor = 5

# Edge mesh spacing discontinuity scaled interpolant and farfield meshing BC
aflr4.input.Mesh_Sizing = {"OML": {"edgeWeight":1.0},
                           "Farfield": {"bcType":"farfield"}}

                           #add scale factor to OML

#------------------------------------------------------------------------------#

# Create AFLR3 AIM to generate the volume mesh
aflr3 = myProblem.analysis.create(aim  = "aflr3AIM",
                                  name = "aflr3")

# Link the aflr4 Surface_Mesh as input to aflr3
aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])


#------------------------------------------------------------------------------#

# Create Fun3D AIM
fun3d = myProblem.analysis.create(aim  = "fun3dAIM",
                                  name = "fun3d")

# Link the aflr3 Volume_Mesh as input to fun3d
fun3d.input["Mesh"].link(aflr3.output["Volume_Mesh"])

# Set project name. Files written to analysisDir will have this name
projectName = "inviscidWing"
fun3d.input.Proj_Name = projectName

fun3d.input.Alpha = 1.0                    # AoA
fun3d.input.Mach = 0.5                     # Mach number
fun3d.input.Equation_Type = "Compressible" # Equation type
fun3d.input.Num_Iter = 5                   # Number of iterations
fun3d.input.Restart_Read = 'off'           # Do not read restart
fun3d.input.Viscous = "inviscid"           # Inviscid calculation

# Set boundary conditions via capsGroup
inviscidBC = {"bcType" : "Inviscid"}
fun3d.input.Boundary_Condition = {"OML"    : inviscidBC,
                                  "Farfield":"farfield"}

# Use python to add inputs to fun3d.nml file
fun3d.input.Use_Python_NML = True

# Write boundary output variables to the fun3d.nml file directly
fun3dnml = f90nml.Namelist()
fun3dnml['boundary_output_variables'] = f90nml.Namelist()
fun3dnml['boundary_output_variables']['mach'] = True
fun3dnml['boundary_output_variables']['cp'] = True
fun3dnml['boundary_output_variables']['average_velocity'] = True

fun3dnml.write(os.path.join(fun3d.analysisDir,"fun3d.nml"), force=True)

# Run AIM pre-analysis
fun3d.preAnalysis()

####### Run FUN3D #####################
print ("\n==> Running FUN3D......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(fun3d.analysisDir) # Move into test directory

# Run fun3d via system call
os.system("nodet_mpi --animation_freq -1 --write_aero_loads_to_file > Info.out")

os.chdir(currentDirectory) # Move back to top directory
#######################################

# Run AIM post-analysis
fun3d.postAnalysis()

# Get force results
print ("\n==> Total Forces and Moments")
# Get Lift and Drag coefficients
print ("--> Cl = ", fun3d.output.CLtot,
           "Cd = ", fun3d.output.CDtot)

# Get Cmx, Cmy, and Cmz coefficients
print ("--> Cmx = ", fun3d.output.CMXtot,
           "Cmy = ", fun3d.output.CMYtot,
           "Cmz = ", fun3d.output.CMZtot)

# Get Cx, Cy, Cz coefficients
print ("--> Cx = ", fun3d.output.CXtot,
           "Cy = ", fun3d.output.CYtot,
           "Cz = ", fun3d.output.CZtot)
