# V3: I think I was forgetting to load the NEXT island's coastline as well, which is necessary for the bridge!

# V2: New idea: insert a point halfway along each artificial bridge, and later use ONLY those 
# as the coastline boundaries of the boxes - ie draw all cross-shore walls starting at one
# of these points, so walls only start in the middle of bridges.  Actually, for the 
# first an last islands in the blob, also have a point at the far side in the middle.  THEN,
# each island will have a north and south box, and no boxes will span two islands.  Boxes might
# be bigger/smaller than 300km*^2 but I'm guessing this is fine.

# Load the coastlines for the "blob" islands and insert artifical 1-D "bridges" between them,
# to be used for defining boxes

# Wait.... since we want an "north" and "south" side of each island, perhaps, for the blobs,
# I should just treat the upper coastline as its own coastline, and the same for the lower.

import os
from pathlib import Path
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
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

output_file_north = os.path.join(output_dir,"coastline_coords_Mercator_island_1_through_2_combined_north")
output_file_south = os.path.join(output_dir,"coastline_coords_Mercator_island_1_through_2_combined_south")

#---------------------------------------------------------------------
#---------------------------------------------------------------------

northern_coast_endpoints_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/s05_blob_coastline_northern_endpoints_per_island.txt"

# load artifical bridge endpoint coordinates
file = open(northern_coast_endpoints_file,'r')
bridge_point_list = file.read().splitlines()
file.close()
bridge_point_list = [ast.literal_eval(el) for el in bridge_point_list]


num_islands_intersecting = 2
coastline_lon_north = []
coastline_lat_north = []
coastline_lon_south = []
coastline_lat_south = []

# ------------------------------------------------------------
# Mods to control the box positions
# ------------------------------------------------------------
forced_box_north_endpoints = []
forced_box_south_endpoints = []
bridge_mid_points_lon = []
bridge_mid_points_lat = []
# ------------------------------------------------------------
first_point_switch = True
# ------------------------------------------------------------


# First, the north coastline of the "blob" islands
for island_number in range(1,num_islands_intersecting+1):   
    
    # Load current island coastline
    coastline_file_in = os.path.join(output_dir,f"coastline_coords_Mercator_island_number_{island_number}.npz")
    d = np.load(coastline_file_in)
    coastline_lonlat = d["coastline_lonlat"]    

    # Load next island coastline
    if island_number < num_islands_intersecting :
        coastline_file_in_next = os.path.join(output_dir,f"coastline_coords_Mercator_island_number_{island_number+1}.npz")
        d = np.load(coastline_file_in)
        coastline_lonlat_next = d["coastline_lonlat"]    
        
    if first_point_switch:
        forced_box_north_endpoints.append(0)        
        first_point_switch = False

    if bridge_point_list[island_number - 1][0] < 0:
        #coastline_lon_north += list(coastline_lonlat[bridge_point_list[island_number - 1][0]-1:,0])
        #coastline_lat_north += list(coastline_lonlat[bridge_point_list[island_number - 1][0]-1:,1])
        coastline_lon_north += list(coastline_lonlat[bridge_point_list[island_number - 1][0]:,0])
        coastline_lat_north += list(coastline_lonlat[bridge_point_list[island_number - 1][0]:,1])
        coastline_lon_north += list(coastline_lonlat[1:bridge_point_list[island_number - 1][1]+1,0])
        coastline_lat_north += list(coastline_lonlat[1:bridge_point_list[island_number - 1][1]+1,1])
        #coastline_lon_north += list(coastline_lonlat[0:bridge_point_list[island_number - 1][1]+1,0])
        #coastline_lat_north += list(coastline_lonlat[0:bridge_point_list[island_number - 1][1]+1,1])
   
    else:
        coastline_lon_north += list(coastline_lonlat[bridge_point_list[island_number - 1][0]:bridge_point_list[island_number - 1][1]+1,0])
        coastline_lat_north += list(coastline_lonlat[bridge_point_list[island_number - 1][0]:bridge_point_list[island_number - 1][1]+1,1])

    # I don't think we need the "middle bridge point" for Mercator.  I think we need a bridge point between islands if there is "more than one grid cell" between them;
    # I guess this needs to be checked by eye (it was needed in WC15) 

#    # Add the new point in the middle of the bridge
#    if island_number < num_islands_intersecting :
#        coastline_lon_north.append((coastline_lonlat[bridge_point_list[island_number - 1][1],0] + coastline_lonlat_next[bridge_point_list[island_number][0],0])/2.0)
#        coast_lat_north.append((coastline_lonlat[bridge_point_list[island_number - 1][1],1] + coastline_lonlat_next[bridge_point_list[island_number][0],1])/2.0)
#       
#        bridge_mid_points_lon.append(coastline_lon_north[-1])
#        bridge_mid_points_lat.append(coast_lat_north[-1])
#
#    # record the index of the added mid-bridge point, for later indexing
#    forced_box_north_endpoints.append(len(coastline_lon_north)-1)
        
coastline_lonlat_north = np.column_stack([coastline_lon_north,coastline_lat_north])


# ------------------------------------------------------------
first_point_switch = True
# ------------------------------------------------------------

# Now, the south coastline of the "blob" islands
#for island_number in range(num_islands_intersecting+1,1,-1):   
for island_number in range(num_islands_intersecting,0,-1):   

    coastline_file_in = os.path.join(output_dir,f"coastline_coords_Mercator_island_number_{island_number}.npz")
    d = np.load(coastline_file_in)
    coastline_lonlat = d["coastline_lonlat"]    

    if first_point_switch:
        forced_box_south_endpoints.append(0)        
        first_point_switch = False

    if bridge_point_list[island_number - 1][0] < 0:
        #coastline_lon_south += list(coastline_lonlat[bridge_point_list[island_number - 1][1]:bridge_point_list[island_number - 1][0],0])
        #coastline_lat_south += list(coastline_lonlat[bridge_point_list[island_number - 1][1]:bridge_point_list[island_number - 1][0],1])
        coastline_lon_south += list(coastline_lonlat[bridge_point_list[island_number - 1][1]:bridge_point_list[island_number - 1][0]+1,0])
        coastline_lat_south += list(coastline_lonlat[bridge_point_list[island_number - 1][1]:bridge_point_list[island_number - 1][0]+1,1])
    
    else:
        coastline_lon_south += list(coastline_lonlat[bridge_point_list[island_number - 1][1]:,0])
        coastline_lat_south += list(coastline_lonlat[bridge_point_list[island_number - 1][1]:,1])
    
#    # Add the new point in the middle of the bridge
#    if island_number >  1 :
#        coastline_lon_south.append(bridge_mid_points_lon.pop())
#        coast_lat_south.append(bridge_mid_points_lat.pop())
#       
#    # record the index of the added mid-bridge point, for later indexing
#    forced_box_north_endpoints.append(len(coastline_lon_north)-1)


coastline_lonlat_south = np.column_stack([coastline_lon_south,coastline_lat_south])


# Plot to see!
fig, ax = plt.subplots()
ax.pcolormesh(lon_field,lat_field,mask,shading="nearest")
ax.plot(coastline_lonlat_north[:,0],coastline_lonlat_north[:,1],c='red')
ax.plot(coastline_lonlat_south[:,0],coastline_lonlat_south[:,1],c='blue')
plt.show()

d = {}
d["coastline_lonlat_north"] = coastline_lonlat_north
np.savez(output_file_north, **d)

d = {}
d["coastline_lonlat_south"] = coastline_lonlat_south
np.savez(output_file_south, **d)







