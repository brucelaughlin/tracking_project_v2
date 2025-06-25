# V2: Changing a few names and targets, so I don't change original files in aa_islands/

# After plotting the isoline with its corresponding "boundary points", I think that it will be easiest
# to proceed if I rotate the isoline until its first entry lines up with the first boundary point.

# BUG!  The original isoline had the starting point tacked onto the end, so it was a closed curve.
# But this was an artificial imposition; I need to un-do it before rotating, and then do the same
# tacking on for the rotated isoline!!!


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

cwd = os.getcwd()
output_dir = os.path.join(str(cwd),"z_output")
Path(output_dir).mkdir(parents=True, exist_ok=True)

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)

lon_field = d["lon_rho"]
lat_field = d["lat_rho"]
mask = d["mask_rho"]

isoline_bp_text_file_out = os.path.join(cwd,'s09_island_isoline_bounding_points_rotated.txt')
#isoline_bp_text_file_out = os.path.join(output_dir,'s09_island_isoline_bounding_points_rotated.txt')

# delete the output text file before writing, to avoid double-writing
os.remove(isoline_bp_text_file_out)

isoline_bp_file_in = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/s07_island_isoline_bounding_points_unrotated.txt"

#---------------------------------------------------------------------
#---------------------------------------------------------------------


# load artifical bounding endpoint coordinates
file = open(isoline_bp_file_in,'r')
isoline_bp_list = file.read().splitlines()
file.close()
isoline_bp_list = [ast.literal_eval(el) for el in isoline_bp_list]


num_islands = 3

num_last_blob_island = 2

island_dex = 0

for island_dex in range(num_last_blob_island,num_islands+1):   
#for island_dex in range(num_last_blob_island,num_last_blob_island+1):   


    if island_dex == num_last_blob_island:
        isoline_file_in = os.path.join(output_dir,'isodistance_lonlat_coords_Mercator_island_1_through_2_blob.npz')
        output_file = os.path.join(output_dir,'isodistance_lonlat_coords_Mercator_island_1_through_2_blob_rotated.npz')
    else:
        isoline_file_in = os.path.join(output_dir,f'isodistance_lonlat_coords_Mercator_island_number_{island_dex}.npz')
        output_file = os.path.join(output_dir,f'isodistance_lonlat_coords_Mercator_island_number_{island_dex}_rotated.npz')

    # Load the isolines
    d = np.load(isoline_file_in)
    isoline = d["isoline_lonlat"]

    isoline_lon = list(isoline[:,0])    
    isoline_lat = list(isoline[:,1])    


    # FIX THE "TACK ON" BUG!!

    # Try... just popping off the tacked-on (ie last) element?!
    isoline_lon.pop()
    isoline_lat.pop()


    isoline_dex = 0    

    while(isoline_dex != isoline_bp_list[island_dex-num_last_blob_island][0]):
        temp_lon = isoline_lon[0]
        temp_lat = isoline_lat[0]
        isoline_lon.pop(0)
        isoline_lat.pop(0)
        isoline_lon.extend([temp_lon])
        isoline_lat.extend([temp_lat])
        isoline_dex += 1

    # Now tack on the starting point as the last element!
    isoline_lon.append(isoline_lon[0])
    isoline_lat.append(isoline_lat[0])

    isoline = np.column_stack([isoline_lon,isoline_lat])


    isoline_bounding_indices = [0,isoline_bp_list[island_dex-num_last_blob_island][1]-isoline_dex]


    # ----------------
    # I had been storing negative values... that gets messy later.  So perhaps I can just make them positive, as follows:
    if isoline_bounding_indices[1] < 0:
        isoline_bounding_indices[1] = isoline_bounding_indices[1] = len(isoline_lon) +  isoline_bounding_indices[1]
    # ----------------

    #with open(isoline_bp_text_file_out, 'a+') as out_file:
    with open(isoline_bp_text_file_out, 'a') as out_file:
        out_file.write(str(isoline_bounding_indices)+'\n')


    d = {}
    d["isoline_lonlat"] = isoline
    np.savez(output_file, **d)










