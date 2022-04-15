import pyCAPS
import os
import argparse

csmFile = "naca_OML.csm"
caps = pyCAPS.Problem(problemName="fluid",
                    capsFile=csmFile,
                    outLevel=1)

#load meshing aims
caps.analysis.create(aim = "egadsTessAIM",name="egadsTess")
meshAIM = caps.analysis.create(aim = "tetgenAIM", name = "tetgen")

#set new egads body tesselation parameters
caps.analysis["egadsTess"].input.Tess_Params = [1.0, 0.01, 20.0]
meshAIM.input.Preserve_Surf_Mesh = True
meshAIM.input["Surface_Mesh"].link(caps.analysis["egadsTess"].output["Surface_Mesh"])

#meshAIM.preAnalysis()
meshAIM.geometry.view()
