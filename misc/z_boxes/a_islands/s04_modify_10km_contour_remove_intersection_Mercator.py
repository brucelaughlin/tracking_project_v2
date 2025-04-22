# Not sure if updating, just using Mercator grid now
# For Mercator, blob islands are 1-2

# OK, new plan, per Chris and Ecologists!  For islands with intersecting isolines, insert artificial
# 0-dimensional coastline "extension" between the islands, remove all of the isoline points in the
# intersections of their polygons, and treat as one big piece of land.

# Furthermore, we want to impose lines dividing "inshore" and "offshore" coasts, per the ecologists.

# V1 is half-working, but need to add the "lower piece" of the 2nd island's isoline FIRST,
# then add the 1st island's isoline, then add the "upper Piece" of the 2nd island's isoline
# (all excluding the isoline points in the intersection)

import os
from pathlib import Path
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from skimage import measure
import scipy.interpolate as spint
from shapely.geometry import Polygon

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


#box_dir = base_path + 'practice/bounding_boxes/create_boxes/'
#islands_dir = 'modify_islands/'
#input_dir = box_dir + islands_dir + 'z_output/'

output_file = os.path.join(output_dir,'isodistance_lonlat_coords_Mercator_island_1_through_2_blob')
#---------------------------------------------------------------------
#---------------------------------------------------------------------

# Function to merge 2 lists (lon/lat) into one list of tuples (GeeksForGeeks):

def merge_lists(list1, list2):
    merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))]
    return merged_list



num_islands_intersecting = 2

island_isolines = []
island_coastlines = []
new_isolines = []


for island_dex in range(1,num_islands_intersecting+1):   

    coastline_file_in = os.path.join(output_dir,f"coastline_coords_Mercator_island_number_{island_dex}.npz")
    isoline_file_in = os.path.join(output_dir,f"isodistance_lonlat_coords_Mercator_island_number_{island_dex}.npz")

    # Load the coastlines
    d = np.load(coastline_file_in)
    coastline_lonlat = d["coastline_lonlat"]    

    island_coastlines.append(coastline_lonlat[0])

    # Load the isolines
    d = np.load(isoline_file_in)
    isoline_lonlat = d["isoline_lonlat"]    
    
    island_isolines.append(isoline_lonlat)


# Need to change things - store the "previous isoline" separately, and 
# only look at the "next" one in the loop
previous_isoline = island_isolines[0] 


for island_dex in range(0,num_islands_intersecting-1):

    combined_isoline_lon = []
    combined_isoline_lat = []

    p = []
    for ii in range(len(previous_isoline)):
        p.append((previous_isoline[ii,0],previous_isoline[ii,1]))

    q = []
    for ii in range(len(island_isolines[island_dex+1])):
        q.append((island_isolines[island_dex+1][ii,0],island_isolines[island_dex+1][ii,1]))


    # Compute the (true mathematical) intersection of the polygons
    pp = Polygon(p)
    qq = Polygon(q)
    internal_points = pp.intersection(qq)

    # Add points not in the intersection to the combined_isoline

    p_lon_pre_1,p_lat_pre_1=internal_points.exterior.xy

    p_lon_pre_2 = list(p_lon_pre_1)    
    p_lat_pre_2 = list(p_lat_pre_1)    

    # Round polygon coordinates and intersection coordinates, otherwise precisions don't match and we can't filter...??
    round_param = 3

    p_lon = list(np.around(np.array(p_lon_pre_2),round_param))
    p_lat = list(np.around(np.array(p_lat_pre_2),round_param))

    coord_tuples_intersection = merge_lists(p_lon,p_lat)

    isoline_2_dex = 0
    # Step 1: Add the "lower piece" of the 2nd island's isoline
    for ii in range(len(island_isolines[island_dex+1])):
        if (round(island_isolines[island_dex+1][ii,0],round_param),round(island_isolines[island_dex+1][ii,1],round_param)) in coord_tuples_intersection:
            isoline_2_dex = ii
            break
        else:
            combined_isoline_lon.append(island_isolines[island_dex+1][ii,0])
            combined_isoline_lat.append(island_isolines[island_dex+1][ii,1])

    # Step 2: Add the 1st island's isoline
    for ii in range(len(previous_isoline)):
        if (round(previous_isoline[ii,0],round_param),round(previous_isoline[ii,1],round_param)) not in coord_tuples_intersection:
            combined_isoline_lon.append(previous_isoline[ii,0])
            combined_isoline_lat.append(previous_isoline[ii,1])

    # Step 3: Add the "upper piece" of the 2nd island's isoline
    for ii in range(isoline_2_dex,len(island_isolines[island_dex+1])):
        if (round(island_isolines[island_dex+1][ii,0],round_param),round(island_isolines[island_dex+1][ii,1],round_param)) not in coord_tuples_intersection:
            combined_isoline_lon.append(island_isolines[island_dex+1][ii,0])
            combined_isoline_lat.append(island_isolines[island_dex+1][ii,1])

    previous_isoline = np.column_stack([combined_isoline_lon,combined_isoline_lat])



# Close and store the final isoline
combined_isoline = np.append(previous_isoline,[previous_isoline[0,:]],axis=0)

#fig, ax = plt.subplots()
#ax.pcolormesh(lon_field,lat_field,mask,shading="nearest")
#ax.plot(combined_isoline[:,0],combined_isoline[:,1])
#plt.show()

d = {}
d["isoline_lonlat"] = combined_isoline
np.savez(output_file, **d)









