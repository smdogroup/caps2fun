"""
Sean Engelstad
test the pip install feature
"""

from pycaps2fun import *

m_namelist = name_index_list(name = ["rib","spar","OML"], index=[10, 2, 30])
print(m_namelist)

aluminum = Material.aluminum()
print(aluminum.dict)

m_shell = ShellProperty(caps_group = "thick1", material=aluminum, thickness=0.2)
print(m_shell.dict)

m_DV = ESP_ThicknessDV(name="thick1", caps_group="thick1",thickness=0.1)
print(m_DV.dv_dict)
