
import os

cwd = os.getcwd()
tacs_dir = os.path.join(cwd, "capsStruct", "Scratch", "tacs")

os.chdir(path=tacs_dir)
for file in os.listdir(tacs_dir):
    if ".f5" in file:
        print(file)
        print(f"\tf5tovtk {file}...")
        os.system(f"~/git/tacs/extern/f5tovtk/f5tovtk {file}")

for file in os.listdir(tacs_dir):
    if ".vtk" in file:
        os.system(f"paraview {file}")

os.chdir(cwd)