# Basic - for sharing box lon/lat with collaborators

# v1 - copied from basic plot in new_island...


#plot_title = 'wc15n model domain\n~300km$^{2}$ coastal boxes (aside from island boxes)\n10km offshore distance as outer wall'
#save_image_name = "domain_full.png"

plot_title = 'Pete Raimandi/Patrick cells'

import netCDF4
import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------------------------------------
#---------------------------------------------------------------------

base_path = '/home/blaughli/tracking_project/'

grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
grid_path_in = base_path + grid_directory + grid_file_in
#grid_path_in = grid_file_in
dset = netCDF4.Dataset(grid_path_in, 'r')

points_type_field = 'rho'
points_type_line = 'psi'
lon_field = np.array(dset['lon_{}'.format(points_type_field)])
lat_field = np.array(dset['lat_{}'.format(points_type_field)])
mask = np.array(dset['mask_{}'.format(points_type_field)])
h = np.array(dset['h'])

dset.close


# Get nice plot background going
# (jet colormap for depth, land masked with grey)

h_2 = np.multiply(mask,h)
cmap_custom = plt.colormaps['jet']
cmap_custom.set_under('0.8')

# ---------------------------------------------



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



fig, ax = plt.subplots()
ax.pcolormesh(lon_field,lat_field,h_2,shading="nearest",cmap = cmap_custom, vmin=0.001)
ax.axis('image')

#print(polygon_vertex_list[0][:,0])

for polygon_dex in range(num_polygons):

    ax.plot(polygon_vertex_list[polygon_dex][:,0],polygon_vertex_list[polygon_dex][:,1],c = 'white',linewidth=0.6)

plt.title(plot_title)

plt.show()


