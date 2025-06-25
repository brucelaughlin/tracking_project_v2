# Hacky
individualIslandsSwitch = True

    import pickle
import os
from pathlib import Path
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import custom_functions.polygon_util as polygon_util
import argparse

visited_points = []

cwd = os.getcwd()
output_dir = os.path.join(str(cwd),"z_output")
Path(output_dir).mkdir(parents=True, exist_ok=True)

output_dir = 'z_output/'

bounding_boxes_file_out = os.path.join(output_dir,"bounding_boxes_lonlat_WCR30_singleCoastalCells.p")
model_file = "/data04/cpennell/runs/wcr30_ERA_v1/wcr30_ERA_20000101/out/roms_avg_wcr30_ERA.nc"

dset = netCDF4.Dataset(model_file,'r')

lon_rho = np.array(dset["lon_rho"][:])
lat_rho = np.array(dset["lat_rho"][:])
mask_rho_original = np.array(dset["mask_rho"][:])

lon_psi = np.array(dset["lon_psi"][:])
lat_psi = np.array(dset["lat_psi"][:])

dset.close()

mask_rho_valid_cells_bool, mask_rho_bool_land_only, mask_rho_valid_cells_keep = polygon_util.generate_walkable_mask(mask_rho_original)

list_of_lists_of_polygon_vertex_coordinates = []

grid_dimensions = np.shape(lon_rho)

# The starting points need to be ocean points JUST OFF the coast - ie start in the first polygon/cell.
# For islands, I start at the first offshore grid cell at the northwest corner of the island.
# note: In WCR30, the three northwest islands in the SCB end up having overlapping coastal polygons
if individualIslandsSwitch:
    starting_point_array = np.array([(1,547),(111,408),(110,414),(112,423),(95,463),(83,463)]) 
else:
    starting_point_array = np.array([(1,547),(111,408),(95,463),(83,463)])

# Update - start from below for continent, from above for island nw corner starting points
starting_point_artificial_previous_point_array = np.copy(starting_point_array)
starting_point_artificial_previous_point_array[0,0] -= 1
starting_point_artificial_previous_point_array[1:,0] += 1

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



