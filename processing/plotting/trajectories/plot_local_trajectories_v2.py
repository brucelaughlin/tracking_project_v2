pdf_file_name = '/home/blaughli/tracking_project/practice/bounding_boxes/final_locations/z_output/binned_data_seasonal_allReleases_baseYear_1999_one_file_pld_45_49_pdrake.npz'
#pdf_file_name = '/home/blaughli/tracking_project/practice/bounding_boxes/final_locations/z_output/z_pre_swap/z_processed_originals/binned_data_seasonal_allReleases_baseYear_1999_one_file_pld_45_49_pdrake.npz'

# -----------------------------------------------------------------------------------------
#save_plot_directory = "/home/blaughli/tracking_project/figures/meetings/will_240702/"
# -----------------------------------------------------------------------------------------

plot_title = 'wc15n model domain\n300km$^{2}$ coastal boxes\n10km offshore distance as outer wall'
save_image_name = "domain_full.png"
# -----------------------------------------------------------------------------------------
#save_plot_file = save_plot_directory + save_image_name
# -----------------------------------------------------------------------------------------



import netCDF4
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
from geopy.distance import great_circle
import scipy.interpolate as spint
import ast
import time

from pyproj import Geod
#from shapely.geometry import Polygon
from shapely.geometry import LineString, Point, Polygon


#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

base_path = '/home/blaughli/tracking_project/'

grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
grid_path_in = base_path + grid_directory + grid_file_in
dset = netCDF4.Dataset(grid_path_in, 'r')

points_type_field = 'rho'
points_type_line = 'psi'
lon_field = np.array(dset['lon_{}'.format(points_type_field)])
lat_field = np.array(dset['lat_{}'.format(points_type_field)])
mask = np.array(dset['mask_{}'.format(points_type_field)])
h = np.array(dset['h'])

dset.close

polygon_file_path = '/home/blaughli/tracking_project/practice/bounding_boxes/final_locations/w_pdrake/s_support_files/wc15.0_06.0km_036km2_settling_polygons.txt'

cell_number = 0
polygon_vertex_list = []
with open(polygon_file_path) as polygon_file:
   for line in polygon_file:
        line_items = line.rstrip().split(',')
        if line_items[0].isdigit():
            if int(line_items[0]) != cell_number:
                if cell_number > 0:
                    polygon_vertex_list.append(current_polygon_vertices)
                cell_number += 1
                current_polygon_vertices = np.array([float(line_items[3]), float(line_items[2])])
                continue
            current_polygon_vertices = np.vstack([current_polygon_vertices, [float(line_items[3]), float(line_items[2])]]) # note that Patrick stores lat first, then lon, so I switch these
# Must append the last polygon
polygon_vertex_list.append(current_polygon_vertices)
num_polygons = len(polygon_vertex_list)



pdf_directory = 'practice/bounding_boxes/final_locations/z_output/z_pre_swap/z_swapped/'

pdf_modified_file = pdf_file_name


d = np.load(pdf_modified_file)

#tick_positions = d['tick_positions']
#tick_labels = d['tick_labels']
#box_num_mod = d['box_num_mod']
local_trajectory_lons = d['local_trajectory_lons']
local_trajectory_lats = d['local_trajectory_lats'] 
local_trajectory_indices = d['local_trajectory_indices']
settle_time_indices = d['settle_time_indices']

pld_days = d['pld_days']

#---------------------------------------------------------------------
#---------------------------------------------------------------------


# Get nice plot background going
# (jet color for depth, land masked with grey)

h_2 = np.multiply(mask,h)
cmap_custom = plt.colormaps['jet']
cmap_custom.set_under('0.8')


# ---------------------------------------------


# I really don't know what "vmin" is doing here, but it seems to work (copied from a stack overflow post)

fig, ax = plt.subplots()
#fig = plt.figure()

ax.pcolormesh(lon_field,lat_field,h_2,shading="nearest",cmap = cmap_custom, vmin=0.001)
ax.axis('image')

for polygon_dex in range(num_polygons):
    ax.plot(polygon_vertex_list[polygon_dex][0],polygon_vertex_list[polygon_dex][1],c = 'white',linewidth=0.6)


# Plot local trajectories for a single box

# Choose a single box...
#chosen_box = 50
chosen_box = 29

num_local = 0

for ii in range(len(local_trajectory_indices)):
    if local_trajectory_indices[ii] == chosen_box:
       ##num_local += 1
        if not np.all(local_trajectory_lons[0:settle_time_indices[ii]] == local_trajectory_lons[ii,0]):
            ax.plot(local_trajectory_lons[ii,0:settle_time_indices[ii]+1],local_trajectory_lats[ii,0:settle_time_indices[ii]+1], c = 'c',linewidth = 1)
            #ax.plot(local_trajectory_lons[ii,0:settle_time_indices[ii]],local_trajectory_lats[ii,0:settle_time_indices[ii]], c = 'c',linewidth = 1)
        #ax.plot(local_trajectory_lons[ii,0:settle_time_indices[ii]],local_trajectory_lats[ii,0:settle_time_indices[ii]], c = 'c',linewidth = 1)
        ax.scatter(local_trajectory_lons[ii,0],local_trajectory_lats[ii,0], c = 'r', s = 10)
        ax.scatter(local_trajectory_lons[ii,settle_time_indices[ii]],local_trajectory_lats[ii,settle_time_indices[ii]], c = 'y', s = 4)


plt.title(plot_title)

#plt.savefig(save_plot_file)
#plt.savefig(save_plot_file, bbox_inches='tight')

#plt.show(bbox_inches='tight')
plt.show()


