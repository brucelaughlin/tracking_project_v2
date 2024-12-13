# Want nice version for reports

# V5: changing the arrows for islands

# V4: just cleaning up; too many commented-out test lines

# V3: remove or improve island labeling for full domain plot

# lat/lon area calc taken from:
# https://stackoverflow.com/questions/68118907/shapely-pyproj-find-area-in-m2-of-a-polygon-created-from-latitude-and-longi


# Using the island coastlines, including the "blob" of islands 1-4, along with the rotated isolines (ie starting points of isolines
# correspond with starting points of coastlines).
# Idea: split each island into two halves, an upper and lower, using the "bounding points" previously determined.  Then proceed
# as if doing the box calculations for two separate coastlines and isolines



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


box_dir = base_path + 'practice/bounding_boxes/create_boxes/'
islands_dir = 'aa_islands/'
continent_dir = 'continent/'
input_dir_islands = box_dir + islands_dir + 'z_output/'
input_dir_continent = box_dir + continent_dir + 'z_output/'

#---------------------------------------------------------------------
#---------------------------------------------------------------------


# Get nice plot background going
# (jet color for depth, land masked with grey)

h_2 = np.multiply(mask,h)
cmap_custom = plt.colormaps['jet']
cmap_custom.set_under('0.8')

# needed?
#plt.ion()


# Now I know the boxes of interest (taken from PDF plots)
tick_positions_orig = [1,5,9,10,13,16,17,19,23,29,32,37,41,47,56,60,67,77]
#tick_positions_orig = [1,5,9,10,13,16,17,19,23,29,32,37,41,47,56,60,67,78]
tick_labels = ['SCl','C','SB','SN','SR','A','SC','SM','TJ','PV','PM','PC','PB','CBl','PR','PA','CM','CB']
tick_positions = [x+1 for x in tick_positions_orig]

island_lats = [32.88,33.36,33.465,33.243,33.967,33.999,34.009,34.039]
island_lons = [-118.483,-118.424,-119.033,-119.518,-120.117,-119.4,-119.774,-120.365]

# ---------------------------------------------
# Plot parameters
box_plot_mod = 2
box_plot_mod_island = 1
#box_plot_mod = 3

font_size_num_continent = 15
font_size_num_island = 12

text_color = "blue"
number_color = "yellow"
#number_color = "white"
arrow_color = "white"
arrow_color2 = "black"
va_val = "center"
offset = 100
#offset = 72
#relpos_tuple = (0,0)
relpos_tuple = (0,0)
relpos_tuple2 = (.5,.5)

# Copied from "https://www.geeksforgeeks.org/matplotlib-pyplot-annotate-in-python/"
bbox = dict(boxstyle ="round", fc ="0.95")
#bbox = dict(boxstyle ="round", fc ="0.8")
arrowprops = dict(
    arrowstyle = "->",
    #relpos = relpos_tuple,
    color=arrow_color,
    connectionstyle = "angle, angleA = 0, angleB = 90, rad = 10")
arrowprops2 = dict(
    arrowstyle = "->",
    #relpos = relpos_tuple,
    color=arrow_color)
arrowprops3 = dict(
    arrowstyle = "->",
    relpos = relpos_tuple,
    color=arrow_color2)
arrowprops4 = dict(
    arrowstyle = "->",
    relpos = relpos_tuple2,
    color=arrow_color2)

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
first_continent_box_dex = 0

# Islands

num_islands = 8
num_last_blob_island = 4

# for labels
box_num = 1
tick_num = 0

# -----------------------------------------------------------------------------------------
# make switch to turn on/off island plotting
switch_plot_islands = True
#switch_plot_islands = False
# -----------------------------------------------------------------------------------------

x_min = -121
x_max = -116.8
y_min = 32.5
y_max = 34.5

# -----------------------------------------------------------------------------------------

for island_dex in range(num_islands,num_last_blob_island-1,-1):
#for island_dex in range(num_last_blob_island,num_islands+1):
#for island_dex in range(num_last_blob_island,num_last_blob_island+1):

    # Set an index for the bounding point lists (the lists of points used to split the coastlines and isolines)
    bp_dex = island_dex-num_last_blob_island

    for inoffshore_switch in range(0,2):


        #if inoffshore_switch == 0:
        if inoffshore_switch == 1:
            bounding_boxes_file_in = input_dir_islands + 'bounding_boxes_lonlat_wc15n_island_number_{}_inshore.p'.format(island_dex)
        else:
            bounding_boxes_file_in = input_dir_islands + 'bounding_boxes_lonlat_wc15n_island_number_{}_offshore.p'.format(island_dex)


        # Load the boxes
        file = open(bounding_boxes_file_in,'rb')
        boxes_lonlat = pickle.load(file)
        file.close

        #for box in boxes_lonlat:
        for box in reversed(boxes_lonlat):
            if box is not None:
                ax.plot(box[0],box[1],c = 'white',linewidth=0.6)
                if box_num == tick_positions[tick_num]:
#tick_positions_orig = [1,5,8,10,13,16,17,19,23,29,32,37,41,47,56,60,67,78]
                    if switch_plot_islands:
                        xy_loc = [island_lons[tick_num], island_lats[tick_num]]
                        #xy_loc = [np.mean(box[0]), np.mean(box[1])]
                        if box_num in [2,11,14]:
                            #ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = [np.mean(box[0]), np.mean(box[1])],
                            ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = xy_loc,
                                xytext =(-.6 * offset, -.2 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops, color=text_color, va=va_val)
                        elif box_num in [6]:
                            ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = xy_loc,
                                xytext =(.4 * offset, -.25 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops, color=text_color, va=va_val)
                        elif box_num in [17]:
                            ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = xy_loc,
                                xytext =(-.5 * offset, -.4 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops, color=text_color, va=va_val)
                        else:
                            ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = xy_loc,
                                xytext =(-.6 * offset, .2 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops, color=text_color, va=va_val)

                    tick_num += 1
                    if tick_num == len(tick_labels):
                        tick_num = 0

                if switch_plot_islands:
                    if box_num % box_plot_mod == 0:
                    #if box_num % box_plot_mod_island == 0:
                        ax.annotate(box_num, xy = [np.mean(box[0]), np.mean(box[1])], color=number_color, ha="center", va="center", fontsize=font_size_num_continent)
                        #ax.annotate(box_num, xy = [np.mean(box[0]), np.mean(box[1])], color=number_color, ha="center", va="center")
                box_num += 1
                first_continent_box_dex += 1
#        plt.show()
#plt.show()

#print("first continent box dex: {}".format(first_continent_box_dex))

#plot_title = '300km^2 coastal boxes\n10km offshore distance as outer wall'
#plt.title(plot_title)
#
#ax.axis('image')
#plt.show()


# Continent

bounding_boxes_file_in = input_dir_continent + 'bounding_boxes_lonlat_coords_{}_coastline_wc15n_continent.p'.format(points_type_line)

# Load the boxes
file = open(bounding_boxes_file_in,'rb')
boxes_lonlat = pickle.load(file)
file.close

for box in boxes_lonlat:
    if box is not None:
        #plt.plot(box[0],box[1],c = 'white',linewidth=0.6)
        ax.plot(box[0],box[1],c = 'white',linewidth=0.6)
        if tick_num < len(tick_labels):
            if box_num == tick_positions[tick_num]:
                #if box_num in [33,38]:
                #if box_num in [30,33,38]:
                if not switch_plot_islands:
                    if box_num in [24]:
                        ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = [np.mean(box[0]), np.mean(box[1])],
                            #xytext =(0, .3 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops2, color=text_color)
                            xytext =(0, .4 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops4, color=text_color)
                    elif box_num in [30]:
                        ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = [np.mean(box[0]), np.mean(box[1])],
                            xytext =(0, .3 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops4, color=text_color)
                    elif box_num in [33]:
                        ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = [np.mean(box[0]), np.mean(box[1])],
                            #xytext =(0, .4 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops3, color=text_color)
                            xytext =(0, .4 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops4, color=text_color)
                    else:
                        ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = [np.mean(box[0]), np.mean(box[1])],
                            #xytext =(.3 * offset, 0), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops2, color=text_color, va=va_val)
                            xytext =(-.7 * offset, 0), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops2, color=text_color, va=va_val)
                tick_num += 1
                if tick_num == len(tick_labels):
                    tick_num = 0
                
            if not switch_plot_islands:
                if (box_num+1) % box_plot_mod == 0:
                #if box_num % box_plot_mod == 0:
                    ax.annotate(box_num, xy = [np.mean(box[0]), np.mean(box[1])], color=number_color, ha="center", va="center", fontsize=font_size_num_continent)
        box_num += 1
        #plt.pause(1)

if switch_plot_islands:
    plt.axis([x_min, x_max, y_min, y_max])


if switch_plot_islands:
    plot_title = '300km^2 coastal boxes\n10km offshore distance as outer wall\nIslands detail'
else:
    plot_title = '300km^2 coastal boxes\n10km offshore distance as outer wall'
plt.title(plot_title)

plt.show()


