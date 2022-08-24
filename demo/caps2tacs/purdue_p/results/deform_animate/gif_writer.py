import utils
import os

gif_writer = utils.GifWriter(frames_per_second=1)
gif_writer(gif_filename="animated.gif", path=os.getcwd())