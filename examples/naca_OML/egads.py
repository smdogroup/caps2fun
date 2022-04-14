import os
import pyCAPS

caps = pyCAPS.Problem(problemName = "struct",
                    capsFile = "naca_OML.csm",
                    outLevel = 1)

wing = caps.geometry

#save cfd egads version
wing.cfgpmtr["cfdOn"].value = 1
wing.save("cfd.egads")

#print("cfdOn Value: ",wing.cfgpmtr["cfdOn"].value)

#save struct egads version
#wing.cfgpmtr["cfdOn"].value = 0
#wing.save("struct.egads")