# Plotting polygons.  For Mercator data, I don't have a grid file, so I'm using Temperature (or maybe another field) as a land-mask substitude.  Hence "crude"

#plot_title = 'wc15n model domain\n~300km$^{2}$ coastal boxes (aside from island boxes)\n10km offshore distance as outer wall'
#save_image_name = "domain_full.png"

plot_title = 'Release/Settlement cells'
#plot_title = 'Pete Raimandi/Patrick cells'

import netCDF4
import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------------------------------------
#---------------------------------------------------------------------

# Mercator file
forcing_file_path = '/data04/cedwards/forcing/mercator/reanalysis12/global-reanalysis-phy-001-030-daily_1993.nc'

#forcing_file_path = grid_file_in
dset = netCDF4.Dataset(forcing_file_path, 'r')

# For Mercator data, lon/lat are 1D lists
lon_field = np.array(dset['longitude'])
lat_field = np.array(dset['latitude'])
T = np.array(dset['thetao'][0,0,:,:])

dset.close

# ---------------------------------------------



polygon_file_path = '/home/blaughli/tracking_project_v2/input_files/wc15.0_06.0km_036km2_settling_polygons.txt'
#polygon_file_path = '/home/blaughli/tracking_project_v2/misc/wc15n_300km2_settling_polygons.txt'

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
ax.pcolormesh(lon_field,lat_field,T)
#ax.pcolormesh(lon_field,lat_field,T,shading="nearest",cmap = cmap_custom, vmin=0.001)
ax.axis('image')

#print(polygon_vertex_list[0][:,0])

for polygon_dex in range(num_polygons):

    ax.plot(polygon_vertex_list[polygon_dex][:,0],polygon_vertex_list[polygon_dex][:,1],c = 'white',linewidth=0.6)

    xy_loc = [np.mean(polygon_vertex_list[polygon_dex][:,0]), np.mean(polygon_vertex_list[polygon_dex][:,1])]

    #ax.annotate(polygon_dex, xy = xy_loc, color='w')


plt.title(plot_title)

plt.show()


