#!/usr/bin/env python3

# Written by Bruce Laughlin, blaughli@ucsc.edu

import sys
import opendrift
from pathlib import Path
import argparse

# Read in the name of the file to process from stdin
parser = argparse.ArgumentParser()
parser.add_argument("trackingfile",type=str)
parser.add_argument("animate",nargs='?',type=str)

args = parser.parse_args()

animate = False

#tracking_output_file = sys.argv[1]
tracking_output_file = args.trackingfile
if args.animate != None:
    animate = True
    print("Creating movie as well as figure")

# Create the output directory "figures" if it doesn't exist already
run_dir = tracking_output_file.rsplit('/',1)[0]
figures_directory = run_dir + '/figures/'
#figures_directory = current_directory + '/figures/'
Path(figures_directory).mkdir(parents=True, exist_ok=True)

# Prepare the names of the output files
#output_file_split = tracking_output_file.split('.')
#output_file_pre = figures_directory + output_file_split[0]
output_file_pre = figures_directory + tracking_output_file.rsplit('/',1)[-1].split('.')[0]
output_png_file = output_file_pre + '.png'
output_mp4_file = output_file_pre + '.mp4'

# Load the model; required for plotting (does magic)
o = opendrift.open(tracking_output_file)

# Create the plot and animation
#o.plot(filename=output_png_file,linecolor="sea_floor_depth_below_sea_level")
o.plot(filename=output_png_file,linecolor="sea_floor_depth_below_sea_level",fast = True)
#o.plot(filename=output_png_file,linecolor="z",fast = True)

if animate:
    o.animation(filename=output_mp4_file,linecolor="z",fast = True)
