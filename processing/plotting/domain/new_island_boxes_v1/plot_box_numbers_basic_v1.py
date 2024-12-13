# Basic - for sharing box lon/lat with collaborators

# v1 - For Pete Raimondi 24/07/23


# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# make switch to turn on/off island plotting
# -----------------------------------------------------------------------------------------

#switch_plot_islands = True
switch_plot_islands = False

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
if switch_plot_islands:
    plot_title = 'wc15n model domain\nNorth/South island coastal boxes\n10km offshore distance as outer wall\n(Southern California Bight detail)'
    save_image_name = "domain_scb.png"
else:
    plot_title = 'wc15n model domain\n~300km$^{2}$ coastal boxes (aside from island boxes)\n10km offshore distance as outer wall'
    save_image_name = "domain_full.png"
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

import netCDF4
import pickle
import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------------------------------------
#---------------------------------------------------------------------

#base_path = '/home/blaughli/tracking_project/'

#grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
#grid_path_in = base_path + grid_directory + grid_file_in
grid_path_in = grid_file_in
dset = netCDF4.Dataset(grid_path_in, 'r')

points_type_field = 'rho'
points_type_line = 'psi'
lon_field = np.array(dset['lon_{}'.format(points_type_field)])
lat_field = np.array(dset['lat_{}'.format(points_type_field)])
mask = np.array(dset['mask_{}'.format(points_type_field)])
h = np.array(dset['h'])

dset.close


##box_dir = base_path + 'practice/bounding_boxes/create_boxes/'
#islands_dir = 'modify_islands/'
#continent_dir = 'continent/'
#input_dir_islands = box_dir + islands_dir + 'z_output/'
#input_dir_continent = box_dir + continent_dir + 'z_output/'


#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Get nice plot background going
# (jet colormap for depth, land masked with grey)

h_2 = np.multiply(mask,h)
cmap_custom = plt.colormaps['jet']
cmap_custom.set_under('0.8')

# ---------------------------------------------

num_islands = 8
num_last_blob_island = 4

# Plot bounds for island plot
x_min = -121
x_max = -116.8
y_min = 32.5
y_max = 34.5


fig, ax = plt.subplots()
ax.pcolormesh(lon_field,lat_field,h_2,shading="nearest",cmap = cmap_custom, vmin=0.001)
ax.axis('image')

for island_dex in range(num_islands,num_last_blob_island-1,-1):

    bounding_boxes_file_in = 'bounding_boxes_lonlat_wc15n_island_number_{}.p'.format(island_dex)
    #bounding_boxes_file_in = input_dir_islands + 'bounding_boxes_lonlat_wc15n_island_number_{}.p'.format(island_dex)

    # Load the boxes
    file = open(bounding_boxes_file_in,'rb')
    boxes_lonlat = pickle.load(file)
    file.close

    for box in reversed(boxes_lonlat):
        if box is not None:
            ax.plot(box[0],box[1],c = 'white',linewidth=0.6)


# Continent

bounding_boxes_file_in = 'bounding_boxes_lonlat_coords_{}_coastline_wc15n_continent.p'.format(points_type_line)
#bounding_boxes_file_in = input_dir_continent + 'bounding_boxes_lonlat_coords_{}_coastline_wc15n_continent.p'.format(points_type_line)

# Load the boxes
file = open(bounding_boxes_file_in,'rb')
boxes_lonlat = pickle.load(file)
file.close

for box in boxes_lonlat:
    if box is not None:
        ax.plot(box[0],box[1],c = 'white',linewidth=0.6)
            

if switch_plot_islands:
    plt.axis([x_min, x_max, y_min, y_max])


plt.title(plot_title)

plt.show()


