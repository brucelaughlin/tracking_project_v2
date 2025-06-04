
import numpy as np
import matplotlib.pyplot as plt

polygon_file_path = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/bounding_boxes_lonlat_Mercator_singleCoastalCells.txt"
grid_file_plot = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_noModification.npz"

d = np.load(grid_file_plot)
mask_rho_plot = d["mask_rho"]
lon_rho = d["lon_rho"]
lat_rho = d["lat_rho"]

cell_number = 0
list_of_polygon_vertex_lonlat_lists = []
with open(polygon_file_path) as polygon_file:
   for line in polygon_file:
        line_items = line.rstrip().split(',')
        if line_items[0].isdigit():
            if int(line_items[0]) != cell_number:
                if cell_number > 0:
                    list_of_polygon_vertex_lonlat_lists.append(current_polygon_vertices)
                cell_number += 1
                current_polygon_vertices = np.array([float(line_items[3]), float(line_items[2])])
                continue
            current_polygon_vertices = np.vstack([current_polygon_vertices, [float(line_items[3]), float(line_items[2])]]) # Lat comes before Lon in the csv file, since that was the Patrick precedent
# Must append the last polygon
list_of_polygon_vertex_lonlat_lists.append(current_polygon_vertices)
num_polygons = len(list_of_polygon_vertex_lonlat_lists)


fig, ax = plt.subplots()
m = ax.pcolormesh(lon_rho,lat_rho,mask_rho_plot,shading="nearest")

for ii in range(len(list_of_polygon_vertex_lonlat_lists)):
    ax.plot(list_of_polygon_vertex_lonlat_lists[ii][:,0],list_of_polygon_vertex_lonlat_lists[ii][:,1], c='r')
#    ax.annotate(ii, xy = [np.mean(list_of_polygon_vertex_lonlat_lists[ii][:,0]), np.mean(list_of_polygon_vertex_lonlat_lists[ii][:,1])], ha="center", va="center", weight="bold",c='k')

ax.axis('image')

plt.show()

