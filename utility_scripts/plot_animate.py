#!/usr/bin/env python3

# Written by Bruce Laughlin, blaughli@ucsc.edu

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from pathlib import Path
import opendrift
import netCDF4
from pathlib import Path

# Create the output directory "figures" if it doesn't exist already
current_directory = os.getcwd()
figures_directory = current_directory + '/figures/'
Path(figures_directory).mkdir(parents=True, exist_ok=True)

# Read in the name of the file to process from stdin
tracking_output_file = sys.argv[1]

# Prepare the names of the output files
output_file_split = tracking_output_file.split('.')
output_file_pre = figures_directory + output_file_split[0]
output_png_file = output_file_pre + '.png'
output_mp4_file = output_file_pre + '.mp4'

# Load the model; required for plotting (does magic)
o = opendrift.open(tracking_output_file)

# Create the plot and animation
o.plot(filename=output_png_file,linecolor="z",fast = True)
o.animation(filename=output_mp4_file,linecolor="z",fast = True)
