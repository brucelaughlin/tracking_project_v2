#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from pathlib import Path
sys.path.append(os.path.abspath("/home/blaughli/tracking_project_v2/python/opendrift_custom/models"))
from larvaldispersal_track_eco_variables import LarvalDispersal
import opendrift
import netCDF4

tracking_output_file = sys.argv[1]

output_file_split = tracking_output_file.split('.')
output_file_pre = output_file_split[0]

#dset = netCDF4.Dataset(tracking_output_file, 'r')
#h_pre = dset['z']
#dset.close
#h = np.array(h_pre)

o = opendrift.open(tracking_output_file)

output_png_file = output_file_pre + '.png'
output_mp4_file = output_file_pre + '.mp4'

o.plot(filename=output_png_file,linecolor="z",fast = True)



