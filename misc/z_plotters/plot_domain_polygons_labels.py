
# lat/lon area calc taken from:
# https://stackoverflow.com/questions/68118907/shapely-pyproj-find-area-in-m2-of-a-polygon-created-from-latitude-and-longi


# Using the island coastlines, including the "blob" of islands 1-4, along with the rotated isolines (ie starting points of isolines
# correspond with starting points of coastlines).
# Idea: split each island into two halves, an upper and lower, using the "bounding points" previously determined.  Then proceed
# as if doing the box calculations for two separate coastlines and isolines


pdf_input_file = "pdf_data_seasonal_ranges_O2_pH_baseYear_1988_WC15N_1988-2010_nRunsPerNode_15_nSeed_020_dielMigration_v1_pld_20_29_swapped.npz"

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
islands_dir = 'modify_islands/'
#islands_dir = 'aa_islands/'
continent_dir = 'continent/'
input_dir_islands = box_dir + islands_dir + 'z_output/'
input_dir_continent = box_dir + continent_dir + 'z_output/'

pdf_directory = 'practice/bounding_boxes/final_locations/z_output/z_pre_swap/z_swapped/'
pdf_modified_file = base_path + pdf_directory + pdf_input_file

#file = open(pdf_modified_file,'rb')
#pdf_list_exposure_T_source_swapped,pdf_list_of_lists_O2_source_swapped,pdf_list_connectivity_swapped,pdf_list_settleTime_swapped,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels,first_continent_box_dex,oxygen_limit_list = pickle.load(file)
#file.close()

d = np.load(pdf_modified_file)

#pdf_arrays_connectivity = d['pdf_arrays_connectivity']
tick_positions = d['tick_positions']
tick_labels = d['tick_labels']
box_num_mod = d['box_num_mod']
#first_continent_box_dex = d['first_continent_box_num']


#---------------------------------------------------------------------
#---------------------------------------------------------------------


# Get nice plot background going
# (jet color for depth, land masked with grey)

h_2 = np.multiply(mask,h)
cmap_custom = plt.colormaps['jet']
cmap_custom.set_under('0.8')

# needed?
#plt.ion()


# Now I know the boxes of interest (taken from domain plots)
#tick_positions_continent_orig = [23,29,32,37,41,47,56,60,67,77]
#tick_labels = ['SCl','Ca','SB','SN','SM','SR','SC','An','TJ','PV','PM','PC','PB','CB','PR','PA','CM','CBl']
#tick_positions = tick_positions_islands + [x+1 for x in tick_positions_continent_orig]

island_lats = [32.823,33.324,33.457,33.224,34.024,33.923,33.957,33.991]
island_lons = [-118.396,-118.345,-119.033,-119.508,-120.338,-120.138,-119.730,-119.4]

#island_lats = [32.823,33.324,33.457,33.224,33.923,33.991,33.957,34.024]
#island_lons = [-118.396,-118.345,-119.033,-119.508,-120.138,-119.4,-119.730,-120.338]

# ---------------------------------------------
# Plot parameters
box_plot_modulo = 2
box_plot_modulo_island = 1
#box_plot_modulo = 3

font_size_num_continent = 15
font_size_num_island = 12
font_size_labels = 17

text_color = "blue"
#number_color = "gold"
#number_color = "mistyrose" # 7
#number_color = "darkorchid"
number_color = "yellow"
#number_color = "white"
#arrow_color = "white"
#arrow_color2 = "darkorange"
#arrow_color2 = "deeppink" # 8
#arrow_color2 = "fuchsia" # 8
arrow_color2 = "magenta" # 9
#arrow_color2 = "violet" # 7
#arrow_color2 = "mediumorchid"
#arrow_color2 = "plum" # 8
#arrow_color2 = "lightgray"
#arrow_color2 = "black"
arrow_width=4
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
    lw=arrow_width,
    color=arrow_color2)
    #color=arrow_color)
arrowprops2 = dict(
    arrowstyle = "->",
    lw=arrow_width,
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

#print(tick_positions)

for island_dex in range(num_islands,num_last_blob_island-1,-1):
#for island_dex in range(num_last_blob_island,num_islands+1):
#for island_dex in range(num_last_blob_island,num_last_blob_island+1):

    # Set an index for the bounding point lists (the lists of points used to split the coastlines and isolines)
    bp_dex = island_dex-num_last_blob_island

    bounding_boxes_file_in = input_dir_islands + 'bounding_boxes_lonlat_wc15n_island_number_{}.p'.format(island_dex)

    # Load the boxes
    file = open(bounding_boxes_file_in,'rb')
    boxes_lonlat = pickle.load(file)
    file.close

    #for box in boxes_lonlat:
    for box in reversed(boxes_lonlat):
        if box is not None:
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


# Continent

bounding_boxes_file_in = input_dir_continent + 'bounding_boxes_lonlat_coords_{}_coastline_wc15n_continent.p'.format(points_type_line)

# Load the boxes
file = open(bounding_boxes_file_in,'rb')
boxes_lonlat = pickle.load(file)
file.close

for box in boxes_lonlat:
    if box is not None:
        ax.plot(box[0],box[1],c = 'white',linewidth=0.6)
        if tick_num < len(tick_labels):
            if box_num == tick_positions[tick_num]:
                if not switch_plot_islands:
                    xy_loc = [np.mean(box[0]), np.mean(box[1])]
                    if box_num in [20]:
                    #if box_num in [24]:
                        ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = xy_loc,
                            xytext =(0, .6 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops, color=text_color, va="center", ha="center",fontsize=font_size_labels)
                    elif box_num in [26]:
                    #elif box_num in [30]:
                        ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = xy_loc,
                            xytext =(0, .5 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops, color=text_color, va="center", ha="center",fontsize=font_size_labels)
                    elif box_num in [31]:
                    #elif box_num in [33]:
                        ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = xy_loc,
                            xytext =(0, .7 * offset), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops, color=text_color, va="center", ha="center",fontsize=font_size_labels)
                    else:
                        ax.annotate("{}: {}".format(box_num,tick_labels[tick_num]), xy = xy_loc,
                            xytext =(-.9 * offset, .0), textcoords ='offset points', bbox = bbox, arrowprops = arrowprops2, color=text_color, va="center", ha="center",fontsize=font_size_labels)
                tick_num += 1
                if tick_num == len(tick_labels):
                    tick_num = 0
                
            #if not switch_plot_islands:
            if (box_num+1) % box_plot_modulo == 0:
            #if box_num % box_plot_modulo == 0:
                ax.annotate(box_num, xy = [np.mean(box[0]), np.mean(box[1])], color=number_color, ha="center", va="center", fontsize=font_size_num_continent, weight="bold")
                #ax.annotate(box_num, xy = [np.mean(box[0]), np.mean(box[1])], color=number_color, ha="center", va="center", fontsize=font_size_num_continent)
        box_num += 1
        #plt.pause(1)

if switch_plot_islands:
    plt.axis([x_min, x_max, y_min, y_max])


plt.title(plot_title)

#plt.savefig(save_plot_file)
#plt.savefig(save_plot_file, bbox_inches='tight')

#plt.show(bbox_inches='tight')
plt.show()


