import capsWrapper

madeupium = capsWrapper.Isotropic.madeupium()
shell_prop = capsWrapper.ShellProperty(caps_group="plate", material=madeupium, membrane_thickness=0.01)
print(shell_prop.dictionary)