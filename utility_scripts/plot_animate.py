#!/usr/bin/env python3

# Written by Bruce Laughlin, blaughli@ucsc.edu

linecolorStr = "z"
#linecolorStr = "sea_floor_depth_below_sea_level"

import os
import opendrift
from pathlib import Path
import argparse

# Process input
parser = argparse.ArgumentParser()
parser.add_argument("trackingfile",type=str)
parser.add_argument("animate",nargs='?',type=str)
args = parser.parse_args()

tracking_output_file = args.trackingfile

tracking_dir = os.path.dirname(os.path.realpath(tracking_output_file))
tracking_file_stem = Path(tracking_output_file).stem

animate = False
if args.animate != None:
    animate = True
    print("Creating movie as well as figure")

# Create the output directory "figures" if it doesn't exist already
figures_dir = os.path.join(tracking_dir,'figures')
Path(figures_dir).mkdir(parents=True, exist_ok=True)

# Prepare the names of the output files
output_png_file_leaf = tracking_file_stem + '.png'
output_mp4_file_leaf = tracking_file_stem + '.mp4'
output_png_file = os.path.join(figures_dir,output_png_file_leaf)
output_mp4_file = os.path.join(figures_dir,output_mp4_file_leaf)

# Create an Opendrift object to access Opendrift plotting methods
o = opendrift.open(tracking_output_file)

# Create the plot and animation
#o.plot(filename=output_png_file,linecolor=linecolorStr)
o.plot(filename=output_png_file,linecolor=linecolorStr,fast = True)

if animate:
    o.animation(filename=output_mp4_file,fast = True)
