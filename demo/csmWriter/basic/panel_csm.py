
from csmWriter import *

csm = CsmFile(filename="panel.csm")

# main file
length = DesignParameter(name="length", value=2.0)
csm.main = length
width = DesignParameter(name="width", value=0.5)
csm.main = width

pattern = ForLoop(num=3)
pattern.statements = Box(xbase=f"3*{pattern.index}", dx=length, dy=width)
csm.main = pattern

csm.main = Sphere(xend=1, yend=1, zend=1)

csm.write()
