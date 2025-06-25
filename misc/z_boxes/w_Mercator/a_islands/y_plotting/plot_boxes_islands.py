# V2: Copied from the domain plots for the report

import netCDF4
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from pathlib import Path


#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"
d = np.load(grid_file)

lon_field = d["lon_rho"]
lat_field = d["lat_rho"]
mask = d["mask_rho"]


box_dir = "/home/blaughli/tracking_project_v2/misc/z_boxes"
islands_dir = 'a_islands/'
continent_dir = 'a_continent/'
input_dir_islands = os.path.join(box_dir,islands_dir,'z_output')
#input_dir_continent = box_dir + continent_dir + 'z_output/'

# Get nice plot background going
# (jet color for depth, land masked with grey)


fig, ax = plt.subplots()
#fig = plt.figure()

ax.pcolormesh(lon_field,lat_field,mask,shading="nearest")
#ax.pcolormesh(lon_field,lat_field,mask,shading="nearest",cmap = cmap_custom, vmin=0.001)
ax.axis('image')

first_continent_box_dex = 0

# Islands

num_islands = 3
num_last_blob_island = 2

# for labels
box_num = 1
tick_num = 0

# Plot bounds for island plot
x_min = -121
x_max = -116.8
y_min = 32.5
y_max = 34.5

# -----------------------------------------------------------------------------------------

for island_dex in range(num_islands,num_last_blob_island-1,-1):
#for island_dex in range(num_last_blob_island,num_islands+1):
#for island_dex in range(num_last_blob_island,num_last_blob_island+1):

    bounding_boxes_file_in = os.path.join(input_dir_islands,'bounding_boxes_lonlat_Mercator_island_number_{}.p'.format(island_dex))

    # Load the boxes
    file = open(bounding_boxes_file_in,'rb')
    boxes_lonlat = pickle.load(file)
    file.close

    #for box in boxes_lonlat:
    for box in reversed(boxes_lonlat):
        if box is not None:
            #print(box_num)
            #print(box)

            ax.plot(box[0],box[1],c = 'black',linewidth=0.6)
            #ax.plot(box[0],box[1],c = 'white',linewidth=0.6)
            ax.annotate(box_num, xy = [np.mean(box[0]), np.mean(box[1])], ha="center", va="center", weight="bold", c='w')
            #first_continent_box_dex += 1
        else:
            print("box {} is empty".format(box_num))
        box_num += 1



# Continent

#bounding_boxes_file_in = input_dir_continent + 'bounding_boxes_lonlat_coords_{}_coastline_Mercator_continent.p'.format(points_type_line)
#
## Load the boxes
#file = open(bounding_boxes_file_in,'rb')
#boxes_lonlat = pickle.load(file)
#file.close
#
#for box in boxes_lonlat:
#    if box is not None:
#        ax.plot(box[0],box[1],c = 'white',linewidth=0.6)
#        ax.annotate(box_num, xy = [np.mean(box[0]), np.mean(box[1])], ha="center", va="center", weight="bold")
#        box_num += 1
#


plt.show()


