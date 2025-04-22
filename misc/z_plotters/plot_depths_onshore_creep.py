
split_index = 1000

import os
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import argparse
from pathlib import Path

## Get the path to the parent directory of the current script
#script_path = os.path.abspath(__file__)
#script_directory = os.path.dirname(script_path)

grid_file = '/home/blaughli/tracking_project_v2/grid_data/wc15n_grd.nc'

parser = argparse.ArgumentParser()
parser.add_argument('inputfile')
args = parser.parse_args()
input_file = args.inputfile

input_file = os.path.abspath(input_file)

input_dir = os.path.dirname(os.path.realpath(input_file))
#input_file_stem = Path(input_file).stem

# Create the output directory "z_figures" if it doesn't exist already
figures_dir = os.path.join(input_dir,'z_figures')
Path(figures_dir).mkdir(parents=True, exist_ok=True)

dset = netCDF4.Dataset(input_file, 'r')
lon = np.array(dset['lon'])
lat = np.array(dset['lat'])
z = np.array(dset['z'])
h_float = -1 * np.array(dset['sea_floor_depth_below_sea_level'])
dset.close()

plot_h_switch = True
if np.std(lon[:,0:-1]) == 0:
    plot_h_switch = False

# The first half of the floats are seeded at the target depth, the 2nd at half the target depth.
# I've decided to plot the shallower seeds, since the deeper ones were seeded in the lowest
# level, with zero velocity
split_index = int(np.shape(z)[0]/2)

#float_index = split_index
float_index = 0

h_float_plot = np.full(np.shape(z)[1], np.nan)

h_float_start_dex = 0
for ii in range(np.shape(z)[1]):
    if h_float[float_index,ii] >= z[float_index,ii]:
        h_float_start_dex = ii
        break

if h_float_start_dex > 3:
    h_float_start_dex -= 4
elif h_float_start_dex > 2:
    h_float_start_dex -= 3
elif h_float_start_dex > 1:
    h_float_start_dex -= 2
elif h_float_start_dex > 0:
    h_float_start_dex -= 1

h_float_plot[h_float_start_dex:] = h_float[float_index,h_float_start_dex:]
h_float_plot = h_float_plot[0:-2]


fig, ax = plt.subplots()

if plot_h_switch:
    ax.plot(h_float_plot, 'r--', label="bathymetry")
ax.plot(z[float_index,:], label="trajectory depth")

title = input_file.split("/")[-2].split("___")[0].split("test_")[-1].split("_seedDepth")[0]
#title = input_file.split("/")[-2].split("___")[0].split("test_")[-1]

ax.set_xlabel('timestep')
ax.set_ylabel('depth (m)')
ax.set_title(title)

plt.legend()
#plt.legend(loc='center right')


# Prepare the names of the output files
#output_png_file_leaf_pre = input_file.split("/")[-2].split("___")[0].split("test_")[-1] + '.png'
output_png_file_leaf = "depth_creep___" + title
output_png_file = os.path.join(figures_dir,output_png_file_leaf)

#plt.show()
plt.savefig(output_png_file)
