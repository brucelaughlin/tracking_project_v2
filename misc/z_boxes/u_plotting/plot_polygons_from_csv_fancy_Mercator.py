# Add labels for polygons with tick marks in connectivity plots

import numpy as np
import matplotlib.pyplot as plt
import argparse

tick_label_csv_file_path = "/home/blaughli/tracking_project_v2/processing/plotting/pdfs/connectivity/t_text_files/tick_labels_single_cell_Mercator.txt"
polygon_file_path = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/bounding_boxes_lonlat_Mercator_singleCoastalCells.txt"
grid_file_plot = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_noModification.npz"

parser = argparse.ArgumentParser()
parser.add_argument("switchscb",type=int)
args = parser.parse_args()

switch_scb = bool(args.switchscb)

#switch_full_domain = True
#switch_scb = False

# Define the tick labels and positions
tick_positions = []
tick_labels = []
with open(tick_label_csv_file_path) as tick_label_file:
   for line in tick_label_file:
        line_items = line.rstrip().split(',')
        if line_items[0].isdigit():
            tick_positions.append(int(line_items[0]))
            tick_labels.append(line_items[1])

tick_positions = np.array(tick_positions)

ticks_to_ignore_full_domain = [-1,91,161,194,181,392]
ticks_to_ignore_scb = [392, tick_positions[-1]+10]


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


font_size_labels = 5

text_color = "blue"
number_color = "yellow"
arrow_color2 = "magenta" # 9
arrow_width=2
#arrow_width=4
offset = 75
#offset = 100

# Copied from "https://www.geeksforgeeks.org/matplotlib-pyplot-annotate-in-python/"
bbox = dict(boxstyle ="round", fc ="0.95")

arrow_props_full_domain = dict(
    arrowstyle = "->",
    lw=arrow_width,
    color=arrow_color2)

lat_max_full_domain = 52
lat_min_full_domain = 31
lon_max_full_domain = -116
lon_min_full_domain = -135

lon_min_scb = -121.8
lon_max_scb = -116.8
lat_min_scb = 32.5
lat_max_scb = 34.5


fig, ax = plt.subplots()

plot_title_pre = "Glorys12 subdomain"

if switch_scb:
    ticks_to_ignore = ticks_to_ignore_scb
    plot_title = f"{plot_title_pre}\nSouthern California Bight detail of study region"
    #plot_title = "Glorys12 connectivity region,\nSouthern California Bight detail"
#if switch_full_domain:
else:
    ticks_to_ignore = ticks_to_ignore_full_domain
    plot_title = f"{plot_title_pre}\nfull study region"
    #plot_title = "Glorys12 connectivity study region"

m = ax.pcolormesh(lon_rho,lat_rho,mask_rho_plot,shading="nearest")

for ii in range(len(list_of_polygon_vertex_lonlat_lists)):
    if switch_scb:
        ax.plot(list_of_polygon_vertex_lonlat_lists[ii][:,0],list_of_polygon_vertex_lonlat_lists[ii][:,1], c='c')
    if ii in tick_positions and ii not in ticks_to_ignore and ii > ticks_to_ignore[0] and ii <= ticks_to_ignore[-1]:
        xy_loc = [np.mean(list_of_polygon_vertex_lonlat_lists[ii][:,0]), np.mean(list_of_polygon_vertex_lonlat_lists[ii][:,1])]
        (ax.annotate("{}: {}".format(ii,tick_labels[np.where(tick_positions==ii)[0][0]]), xy = xy_loc,
            xytext =(-.9 * offset, .0), textcoords ='offset points', bbox = bbox, 
        arrowprops = arrow_props_full_domain, color=text_color, va="center", ha="center",fontsize=font_size_labels))
    

ax.axis('image')

if switch_scb:
    plt.xlim(lon_min_scb,lon_max_scb)
    plt.ylim(lat_min_scb,lat_max_scb)
#if switch_full_domain:
else:
    plt.xlim(lon_min_full_domain,lon_max_full_domain)
    plt.ylim(lat_min_full_domain,lat_max_full_domain)

plt.title(plot_title)

plt.show()

