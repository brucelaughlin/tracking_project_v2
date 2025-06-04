

tick_label_csv_file_path = "/home/blaughli/tracking_project_v2/processing/plotting/pdfs/connectivity/t_text_files/tick_labels_single_cell_Mercator.txt"


# -----------------------------------------------------------------------------------------
save_plot_directory = "/home/blaughli/tracking_project/figures/meetings/will_240702/"
# -----------------------------------------------------------------------------------------




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
    #plot_title = 'wc15n model domain\n300km$^{2}$ coastal boxes\n10km offshore distance as outer wall\n(Southern California Bight detail)'
    save_image_name = "domain_scb.png"
else:
    plot_title = 'wc15n model domain\n300km$^{2}$ coastal boxes\n10km offshore distance as outer wall'
    save_image_name = "domain_full.png"
# -----------------------------------------------------------------------------------------
save_plot_file = save_plot_directory + save_image_name
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



# Define the tick labels and positions
tick_positions = []
tick_labels = []
with open(tick_label_csv_file_path) as tick_label_file:
   for line in tick_label_file:
        line_items = line.rstrip().split(',')
        if line_items[0].isdigit():
            tick_positions.append(int(line_items[0]))
            tick_labels.append(line_items[1])






font_size_num_continent = 15
font_size_num_island = 12
font_size_labels = 17

text_color = "blue"
number_color = "yellow"
arrow_color = "magenta"
arrow_width=4
va_val = "center"
offset = 100
#offset = 72
relpos_tuple = (0,0)
relpos_tuple2 = (.5,.5)

# Copied from "https://www.geeksforgeeks.org/matplotlib-pyplot-annotate-in-python/"
bbox = dict(boxstyle ="round", fc ="0.95")
#bbox = dict(boxstyle ="round", fc ="0.8")
arrowprops = dict(
    arrowstyle = "->",
    lw=arrow_width,
    color=arrow_color)
arrowprops2 = dict(
    arrowstyle = "->",
    lw=arrow_width,
    color=arrow_color)

# ---------------------------------------------


# I really don't know what "vmin" is doing here, but it seems to work (copied from a stack overflow post)

fig, ax = plt.subplots()
#fig = plt.figure()

ax.pcolormesh(lon_field,lat_field,h_2,shading="nearest",cmap = cmap_custom, vmin=0.001)
#plt.pcolormesh(lon_field,lat_field,h_2,shading="nearest",cmap = cmap_custom, vmin=0.001)
ax.axis('image')
#fig.axis('image')

#fig.canvas.draw()
#plt.show()

# for figuring out the number of the last island box
#first_continent_box_dex = 0

# Islands

num_islands = 8
num_last_blob_island = 4

# for labels
box_num = 1
tick_num = 0

# Plot bounds for island plot
x_min = -121
x_max = -116.8
y_min = 32.5
y_max = 34.5

# -----------------------------------------------------------------------------------------


fig, ax = plt.subplots()
m = ax.pcolormesh(lon_rho,lat_rho,mask_rho_plot,shading="nearest")

for ii in range(len(list_of_polygon_vertex_lonlat_lists)):
    ax.plot(list_of_polygon_vertex_lonlat_lists[ii][:,0],list_of_polygon_vertex_lonlat_lists[ii][:,1], c='r')
#    ax.annotate(ii, xy = [np.mean(list_of_polygon_vertex_lonlat_lists[ii][:,0]), np.mean(list_of_polygon_vertex_lonlat_lists[ii][:,1])], ha="center", va="center", weight="bold",c='k')







#print(tick_positions)

            ax.plot(box[0],box[1],c = 'white',linewidth=0.6)
            #if (box_num-1) == tick_positions[tick_num]:
            if box_num == tick_positions[tick_num]:
#tick_positions_orig = [1,5,9,10,13,16,17,19,23,29,32,37,41,47,56,60,67,77]
                if switch_plot_islands:
                    xy_loc = [island_lons[tick_num], island_lats[tick_num]]

                    #ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = xy_loc,
                    ax.annotate("{}: {}".format(tick_positions[tick_num],tick_labels[tick_num]), xy = xy_loc,
                    #ax.annotate("{}: {}".format(box_numbers_islands_mod[box_num-1]+1,tick_labels[tick_num]), xy = xy_loc,
                        xytext =(0, -.8 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops, color=text_color, va="center", ha="center",fontsize=font_size_labels)

                tick_num += 1
                if tick_num == len(tick_labels):
                    tick_num = 0

            if switch_plot_islands:
                #if box_num % box_plot_modulo == 0:
                if box_num % box_plot_modulo_island == 0:
                    #ax.annotate(box_numbers_islands_mod[box_num-1]+1, xy = [np.mean(box[0]), np.mean(box[1])], color=number_color, ha="center", va="center", fontsize=font_size_num_continent, weight="bold")
                    ax.annotate(box_num_mod[box_num-1]+1, xy = [np.mean(box[0]), np.mean(box[1])], color=number_color, ha="center", va="center", fontsize=font_size_num_continent, weight="bold")
                    #ax.annotate(box_num, xy = [np.mean(box[0]), np.mean(box[1])], color=number_color, ha="center", va="center", fontsize=font_size_num_continent, weight="bold")
            box_num += 1
            #first_continent_box_dex += 1









if switch_plot_islands:
    plt.axis([x_min, x_max, y_min, y_max])


plt.title(plot_title)

#plt.savefig(save_plot_file)
#plt.savefig(save_plot_file, bbox_inches='tight')

#plt.show(bbox_inches='tight')
plt.show()


