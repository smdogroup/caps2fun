import os
import numpy as np
from fileIO import writeInput, makeDVdict, readOutput

#set the functions to check
functions = ["ksfailure","cl","cd","mass"]

#generate variables and functions from writeInput.py
DVdict = makeDVdict()

#settings
#f2f analysis: "aerothermal", "aeroelastic", "aerothermoelastic"
#fun3d analysis: "inviscid", "laminar", "turbulent"
#modes: "forward", "adjoint", "complex_step"

#write the F2F input file for complex mode
writeInput(DVdict, functions, f2fanalysis="aeroelastic", fun3danalysis="laminar",  mode="forward")