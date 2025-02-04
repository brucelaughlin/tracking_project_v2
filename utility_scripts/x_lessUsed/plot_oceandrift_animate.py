#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from pathlib import Path
sys.path.append(os.path.abspath("/home/blaughli/tracking_project/opendrift_custom/models"))
import opendrift
import netCDF4

#base_dir = '/home/blaughli/tracking_project/practice/experiments/practice_1/sanity_check/z_output/'


tracking_output_file = sys.argv[1]


#itracking_output_file = 'test_output_01_01_01_run_1_of_1.nc'
#o = opendrift.open(tracking_output_file)

output_file_split = tracking_output_file.split('.')
#output_file_pre = 'z_media/' + output_file_split[0]
output_file_pre = output_file_split[0]

dset = netCDF4.Dataset(tracking_output_file, 'r')

h_pre = dset['z']
h = np.array(h_pre)

o = opendrift.open(tracking_output_file)

#output_png_file = 'overhead_plot_all.png'
#output_z_png_file = 'z_plot_all.png'
#output_mp4_file = 'overhead_animation_all.mp4'

output_png_file = output_file_pre + '.png'
output_mp4_file = output_file_pre + '.mp4'

#o.plot(filename=output_png_file)
o.plot(filename=output_png_file,linecolor="z",fast = True)
#o.plot_property('z',filename=output_z_png_file)
o.animation(filename=output_mp4_file,linecolor="z",fast = True)



