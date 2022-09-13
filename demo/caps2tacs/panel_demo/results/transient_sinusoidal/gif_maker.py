import os
import utils

# all the capsWrapper package modules/classes are also available 
# in caps2tacs since caps2tacs package imports the whole capsWrapper package\

gif_writer = utils.GifWriter(frames_per_second=20)
gif_writer(gif_filename="transient_sinusoidal.gif", path=os.getcwd())


