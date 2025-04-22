
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


###fig, ax = plt.subplots()

#fig, host = plt.subplots(figsize=(8,5), layout='constrained')
fig, host = plt.subplots()

fig.subplots_adjust(bottom=0.2, left=0.2, right=0.8)

title = "depth and ocean vertical diffusivity along trajectory"
footTitle = input_file.split("/")[-2].split("___")[0].split("test_")[-1].split("_seedDepth")[0]
#title = input_file.split("/")[-2].split("___")[0].split("test_")[-1]



#ax.set_xlabel(f'timestep\n{footTitle}')
#ax.set_ylabel(r'diffusivity ($m^2/s$)')
#ax.set_title(title)

ax2 = host.twinx()

host.set_xlabel(f'timestep\n\n{footTitle}')
host.set_ylabel("depth (m)")
host.set_title(title)

ax2.set_ylabel(r'diffusivity ($m^2/s$)')


p1 = host.plot(z[float_index,:], label="trajectory depth", color='blue')

p2 = ax2.plot(aks[float_index,0:-2], label="vertical diffusivity", color='red')
#p2 = ax2.plot(aks[float_index,0:-2], 'r--', label="vertical diffusivity")

#host.legend(handles=p1+p2, loc='best')

ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

#plt.legend()
###plt.legend(loc='center right')

host.yaxis.label.set_color(p1[0].get_color())
ax2.yaxis.label.set_color(p2[0].get_color())


# Prepare the names of the output files
#output_png_file_leaf_pre = input_file.split("/")[-2].split("___")[0].split("test_")[-1] + '.png'
output_png_file_leaf = "vertical_diffusivity___" + footTitle
output_png_file = os.path.join(figures_dir,output_png_file_leaf)

#plt.show()
plt.savefig(output_png_file)
