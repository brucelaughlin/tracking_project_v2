# Inputs:
# <modelfile> : A single file from which to load grid data; either a ROMS netcdf file or an npz file of my own creation (derived from another model, ie Mercator, etc)
# <startingpointfile> : A newline-separated file whose lines contain single tuples of (i,j) starting points for walking along coasts while storing single-cell polygons.
# Typically I put these starting points at the first Northwest polygon of islands, and at the first Southwest polygon of a western coastline
# <nonromsswitch> : If anything is provided, the model is assumed to not be ROMS, and the grid data will be loaded from a custom-made ".npz" file of my own creation

import pickle
import os
from pathlib import Path
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import custom_functions.polygon_util as polygon_util
import argparse
import ast

parser = argparse.ArgumentParser()
parser.add_argument("modelfile", type=str)
parser.add_argument("startingpointfile", type=str)
parser.add_argument("nonromsswitch", type=str, nargs="?")
args = parser.parse_args()

model_file = args.modelfile
starting_point_file = args.startingpointfile
non_roms_switch = args.nonromsswitch

cwd = os.getcwd()
output_dir = os.path.join(str(cwd),"z_output")
Path(output_dir).mkdir(parents=True, exist_ok=True)

starting_point_list = []

with open(starting_point_file, 'r') as infile:
    for line in infile:
        starting_point_list.append(ast.literal_eval(line))

starting_point_array = np.array(starting_point_list)

model_label = Path(model_file).stem

bounding_boxes_file_out_pre = f"bounding_boxes_lonlat__{model_label}__singleCoastalCells.p"
bounding_boxes_file_out = os.path.join(output_dir, bounding_boxes_file_out_pre)

if non_roms_switch is not None:
    d = np.load(model_file)
    mask_rho_original = d["mask_rho"]
    lon_rho = d["lon_rho"]
    lat_rho = d["lat_rho"]
    lon_psi = d["lon_psi"]
    lat_psi = d["lat_psi"]
else:
    dset = netCDF4.Dataset(model_file,'r')
    mask_rho_original = np.array(dset["mask_rho"][:])
    lon_rho = np.array(dset["lon_rho"][:])
    lat_rho = np.array(dset["lat_rho"][:])
    lon_psi = np.array(dset["lon_psi"][:])
    lat_psi = np.array(dset["lat_psi"][:])
    dset.close()

mask_rho_valid_cells_bool, mask_rho_bool_land_only, mask_rho_valid_cells_keep = polygon_util.generate_walkable_mask(mask_rho_original)

list_of_lists_of_polygon_vertex_coordinates = []

grid_dimensions = np.shape(lon_rho)

# Update - start from below for continent, from above for island nw corner starting points
starting_point_artificial_previous_point_array = np.copy(starting_point_array)
starting_point_artificial_previous_point_array[0,0] -= 1
starting_point_artificial_previous_point_array[1:,0] += 1

visited_points = []

for i_starting_point in range(len(starting_point_array)):

    polygon_util.dfs(tuple(starting_point_array[i_starting_point]),tuple(starting_point_artificial_previous_point_array[i_starting_point]), visited_points, mask_rho_valid_cells_bool, mask_rho_bool_land_only)


polygon_util.generate_single_cell_polygons_from_list_of_points(list_of_lists_of_polygon_vertex_coordinates,visited_points,lon_psi,lat_psi)

fig, ax = plt.subplots()
m = ax.pcolormesh(lon_rho,lat_rho,mask_rho_valid_cells_keep,shading="nearest")

for ii in range(len(list_of_lists_of_polygon_vertex_coordinates)):
    ax.plot(list_of_lists_of_polygon_vertex_coordinates[ii][0,:],list_of_lists_of_polygon_vertex_coordinates[ii][1,:], c='r')
    ax.annotate(ii, xy = [np.mean(list_of_lists_of_polygon_vertex_coordinates[ii][0,:]), np.mean(list_of_lists_of_polygon_vertex_coordinates[ii][1,:])], ha="center", va="center", weight="bold",c='k')

ax.axis('image')

plt.show()


file = open(bounding_boxes_file_out,'wb')
pickle.dump(list_of_lists_of_polygon_vertex_coordinates,file)
file.close()



