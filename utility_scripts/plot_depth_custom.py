#!/usr/bin/env python3

# Stolen from "__init__.py" -> "plot_property"

import sys
import os
from pathlib import Path
import opendrift
import matplotlib.pyplot as plt
from matplotlib import dates
import netCDF4
import numpy as np

# Create the output directory "figures" if it doesn't exist already
current_directory = os.getcwd()
figures_directory = current_directory + '/figures/'
Path(figures_directory).mkdir(parents=True, exist_ok=True)

# Read in the name of the file to process from stdin
tracking_file = sys.argv[1]


output_figures_dir_stem = "z_figures"

tracking_dir = os.path.dirname(os.path.realpath(tracking_file))
tracking_file_stem = Path(tracking_file).stem


# Create the output directories if they don't exist already
#base = os.path.splitext(tracking_file)[0]

output_figures_dir = os.path.join(tracking_dir,output_figures_dir_stem)
#output_figures_dir = base.rsplit('/', 1)[0] + "/" + output_figures_dir_stem
Path(output_figures_dir).mkdir(parents=True, exist_ok=True)


figure_file_leaf = tracking_file_stem + '_depth_custom.png'
figure_file = os.path.join(output_figures_dir,figure_file_leaf)


# Prepare the names of the output files
#output_file_split = tracking_file.split('.')
#output_file_pre = figures_directory + output_file_split[0]
#output_png_file = output_file_pre + '_depth_custom.png'

dset = netCDF4.Dataset(tracking_file)
z = dset.variables['z'][:].data
bottom_depth = dset.variables['sea_floor_depth_below_sea_level'][:].data
dset.close()

initial_bottom_depths = []
for trajectory_dex in range(np.shape(bottom_depth)[0]):
    initial_bottom_depths.append(int(bottom_depth[trajectory_dex,0]))
    #initial_bottom_depths.append(bottom_depth[trajectory_dex,0])


"""Basic function to plot time series of any element properties."""

fig = plt.figure()
ax = fig.gca()
plt.xticks(rotation='vertical')


for trajectory_dex in range(np.shape(bottom_depth)[0]):
    plt.plot(z[trajectory_dex,:], label=f"bottom_depth: {bottom_depth[trajectory_dex,0]}m")

plt.title('z')
plt.ylabel('z (m)')
plt.xlabel('timestep (hour)')
plt.subplots_adjust(bottom=.3)
plt.grid()
plt.legend()

if figure_file is None:
    plt.show()
else:
    plt.savefig(figure_file)

