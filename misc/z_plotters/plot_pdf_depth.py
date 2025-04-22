target_depth = 1000


output_filename_pre = "trajectory_depths_pdf_"

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
aks = np.array(dset['ocean_vertical_diffusivity'])
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

# The last element is a NaN, and the 2nd to last might also be problematic??
z_toBin = z[0,0:-2] + target_depth
#h_float_toBin = h_float[0,0:-2]


#hist, bin_edges = np.histogram(h_float_plot, bins=10, density=True)


###fig, ax = plt.subplots()

#fig, ax = plt.subplots(figsize=(8,5), layout='constrained')
fig, ax = plt.subplots()

fig.subplots_adjust(bottom=0.2, left=0.2)
#fig.subplots_adjust(bottom=0.2, left=0.2, right=0.8)

title = "pdf of depth under mixing, behavior, or both"
footTitle = input_file.split("/")[-2].split("___")[0].split("test_")[-1].split("_seedDepth")[0]
#title = input_file.split("/")[-2].split("___")[0].split("test_")[-1]


ax.set_xlabel(f'depth relative to target depth (m)\n\n{footTitle}')
ax.set_ylabel("probability")
ax.set_title(title)

#p1 = ax.plot(z[float_index,:], label="trajectory depth", color='blue')

#ax.hist(z_toBin)
#ax.hist(z_toBin, density=True)
ax.hist(z_toBin, bins=25, density=True)

# Prepare the names of the output files
#output_png_file_leaf_pre = input_file.split("/")[-2].split("___")[0].split("test_")[-1] + '.png'
output_png_file_leaf = output_filename_pre + footTitle
output_png_file = os.path.join(figures_dir,output_png_file_leaf)

#plt.show()
plt.savefig(output_png_file)
