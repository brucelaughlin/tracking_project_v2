# V2: I think the last point in the north coast and the first point in the south coast are the same.  For my
# new box algorithm, with fixed cross-shore locations, I think it'll be easier ot just have one coastline file.
# So, I want to remove the redundant point and then merge the north and south, to create one coastline list.

import os
from pathlib import Path
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from skimage import measure
import scipy.interpolate as spint
import ast


#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#cwd = os.getcwd()
#output_dir = os.path.join(str(cwd),"z_output")
#Path(output_dir).mkdir(parents=True, exist_ok=True)
output_dir = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/z_output"

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)

lon_field = d["lon_rho"]
lat_field = d["lat_rho"]
mask = d["mask_rho"]

#---------------------------------------------------------------------
#---------------------------------------------------------------------

# load artifical bounding endpoint coordinates

island_bounding_point_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/y_plotting/v_final_island_coastline_bounding_points.txt"

file = open(island_bounding_point_file,'r')
bounding_point_list = file.read().splitlines()
file.close()
bounding_point_list = [ast.literal_eval(el) for el in bounding_point_list]



fig, ax = plt.subplots()
ax.pcolormesh(lon_field,lat_field,mask,shading="nearest")

num_islands = 3

num_last_blob_island = 2

island_dex = 0


coastline_inshore_file_in = os.path.join(output_dir,'coastline_coords_Mercator_island_1_through_2_combined_north.npz')
coastline_offshore_file_in = os.path.join(output_dir,'coastline_coords_Mercator_island_1_through_2_combined_south.npz')

isoline_file_in = os.path.join(output_dir,'isodistance_lonlat_coords_Mercator_island_1_through_2_blob_rotated.npz')
#isoline_file_in = os.path.join(output_dir,'isodistance_lonlat_coords_coastline_Mercator_island_1_through_2_blob_rotated.npz')

# Load the coastlines
d = np.load(coastline_inshore_file_in)
coastline_lonlat_1 = d["coastline_lonlat_north"]

d = np.load(coastline_offshore_file_in)
coastline_lonlat_2 = d["coastline_lonlat_south"]
coastline = np.append(coastline_lonlat_1,coastline_lonlat_2[1:,:],axis=0)
