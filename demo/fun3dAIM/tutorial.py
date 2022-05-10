#demo for fun3d aim
# Create Fun3D AIM
fun3d = myProblem.analysis.create(aim = "fun3dAIM",
name = "fun3d")
# Link the aflr3 Volume_Mesh as input to fun3d
fun3d.input["Mesh"].link(aflr3.output["Volume_Mesh"])
# Set project name. Files written to analysisDir will have this name
projectName = "inviscidWing"
fun3d.input.Proj_Name = projectName
fun3d.input.Alpha = 1.0 # AoA
fun3d.input.Mach = 0.5 # Mach number
fun3d.input.Equation_Type = "Compressible" # Equation type
fun3d.input.Num_Iter = 5 # Number of iterations
fun3d.input.Restart_Read = ’off’ # Do not read restart
fun3d.input.Viscous = "inviscid" # Inviscid calculation
# Set boundary conditions via capsGroup
inviscidBC = {"bcType" : "Inviscid"}
fun3d.input.Boundary_Condition = {"Wing" : inviscidBC,
"Farfield":"farfield"}