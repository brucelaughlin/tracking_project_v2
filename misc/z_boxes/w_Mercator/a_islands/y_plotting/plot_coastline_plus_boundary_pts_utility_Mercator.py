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
        
for island_dex in range(num_last_blob_island,num_islands+1):   
#for island_dex in range(num_last_blob_island,num_last_blob_island+1):   


    if island_dex == num_last_blob_island:
    
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

    else:
        coastline_file_in = os.path.join(output_dir,'coastline_coords_Mercator_island_number_{}.npz').format(island_dex)
        isoline_file_in = os.path.join(output_dir,'isodistance_lonlat_coords_Mercator_island_number_{}_rotated.npz').format(island_dex)
        #isoline_file_in = os.path.join(output_dir,'isodistance_lonlat_coords_coastline_Mercator_island_number_{}_rotated.npz').format(island_dex)

        # Load the coastlines
        d = np.load(coastline_file_in)
        coastline_lonlat = d["coastline_lonlat"]
        coastline = coastline_lonlat
        #coastline = coastline_lonlat[0]

    # Load the isolines
    d = np.load(isoline_file_in)
    isoline_lonlat = d["isoline_lonlat"]

    isoline = isoline_lonlat

    ax.plot(coastline[:,0],coastline[:,1])
    ax.plot(isoline[:,0],isoline[:,1])
    ax.scatter(coastline[0,0],coastline[0,1],c='blue')

    ax.scatter(coastline[bounding_point_list[island_dex-num_last_blob_island][0],0],coastline[bounding_point_list[island_dex-num_last_blob_island][0],1],c='red')
    ax.scatter(coastline[bounding_point_list[island_dex-num_last_blob_island][1],0],coastline[bounding_point_list[island_dex-num_last_blob_island][1],1],c='red')


    for ii in range(np.shape(coastline)[0]):
        ax.annotate(ii, xy = [coastline[ii,0],coastline[ii,1]], color='white', ha="center", va="center", fontsize=15, weight="bold")
    for ii in range(np.shape(isoline)[0]):
        ax.annotate(ii, xy = [isoline[ii,0],isoline[ii,1]], color='white', ha="center", va="center", fontsize=15, weight="bold")



    #ax.axis('image')
    #plt.show()
ax.axis('image')
plt.show()









