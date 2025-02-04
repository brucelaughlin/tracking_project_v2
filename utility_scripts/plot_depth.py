#!/usr/bin/env python3

import sys
import os
from pathlib import Path
import opendrift

# Create the output directory "figures" if it doesn't exist already
current_directory = os.getcwd()
figures_directory = current_directory + '/figures/'
Path(figures_directory).mkdir(parents=True, exist_ok=True)

# Read in the name of the file to process from stdin
tracking_output_file = sys.argv[1]

# Prepare the names of the output files
output_file_split = tracking_output_file.split('.')
output_file_pre = figures_directory + output_file_split[0]
output_png_file = output_file_pre + '_depth.png'

# Load the model; required for plotting (does magic)
o = opendrift.open(tracking_output_file)

# Create the plot
#o.plot_property('z',filename=output_png_file,fast = True)
o.plot_property('z',filename=output_png_file)



