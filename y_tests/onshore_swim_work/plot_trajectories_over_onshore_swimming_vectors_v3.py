# v1: trying to figure out why particles are being deactivated, specifically in the gradient swimming runs.  So, plotting
# over the swim mask, or can plot over rho mask

# https://stackoverflow.com/questions/20165169/change-colour-of-curve-according-to-its-y-value-in-matplotlib

# Input Files
#---------------------------------------------------------------------

#tracking_output_dir_pre = "drift_150_swim_test_v2_zeroVel_grad/"
#swim_data_file = '/home/blaughli/tracking_project/practice/onshore_swim_work/z_production/swim_data.p'

tracking_output_dir_pre = "drift_150_swim_test_v2_zeroVel/"
swim_data_file = '/home/blaughli/tracking_project/practice/onshore_swim_work/z_production/swim_data_simple_v2.npz'

#---------------------------------------------------------------------




#particle_numbers_islands = list(range(10000,50000,10000))
particle_numbers_full = list(range(78))

#colors_paths = ['skyblue','lime','violet','orange']


# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# make switch to turn on/off island plotting
# -----------------------------------------------------------------------------------------
#switch_plot_islands = True
switch_plot_islands = False
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
#if switch_plot_islands:
#    plot_title = 'wc_15n model domain\n300km$^{2}$ coastal boxes\n10km offshore distance as outer wall\n(Southern California Bight detail)'
#    save_image_name = "domain_scb.png"
#else:
#    plot_title = 'wc_15n model domain\n300km$^{2}$ coastal boxes\n10km offshore distance as outer wall'
#    save_image_name = "domain_full.png"


fig_paramTitle = "wc_15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"
fig_mainTitle = "Sample trajectories"

plot_title = fig_mainTitle + "\n" + fig_paramTitle

tracking_output_file = '/data/blaughli/tracking_output/swim_tests/dummy_dir_onshoreSwim_test/tracking_output_calcDT_060_saveDT_1440_buffer_100_nSeed_001_startNudge_000000.nc'

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
from os import listdir
from os.path import isfile, join
import sys
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.lines import Line2D

from pyproj import Geod
#from shapely.geometry import Polygon
from shapely.geometry import LineString, Point, Polygon



#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#tracking_output_base = "/data01/blaughli/tracking_project_output/"
#tracking_output_dir = tracking_output_base + tracking_output_dir_pre

base_path = '/home/blaughli/tracking_project/'

grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
grid_path_in = base_path + grid_directory + grid_file_in
dset = netCDF4.Dataset(grid_path_in, 'r')

points_type_field = 'rho'
points_type_line = 'psi'
lon_field = np.array(dset['lon_{}'.format(points_type_field)])
lat_field = np.array(dset['lat_{}'.format(points_type_field)])
mask_field = np.array(dset['mask_{}'.format(points_type_field)])
h = np.array(dset['h'])

dset.close


box_dir = base_path + 'practice/bounding_boxes/create_boxes/'
islands_dir = 'modify_islands/'
continent_dir = 'continent/'
input_dir_islands = box_dir + islands_dir + 'z_output/'
input_dir_continent = box_dir + continent_dir + 'z_output/'


#tracking_output_files = [f for f in listdir(tracking_output_dir) if isfile(join(tracking_output_dir,f))]
#tracking_output_files.sort()


d = np.load(swim_data_file)
onshore_swim_component_x_map = d['onshore_swim_component_x_map']
onshore_swim_component_y_map = d['onshore_swim_component_y_map']
mask = d['mask']

#file = open(swim_data_file,'rb')
#mask,mask_flat,coord_array,onshore_swim_component_x,onshore_swim_component_y,onshore_swim_component_x_map,onshore_swim_component_y_map = pickle.load(file)
#file.close()


#---------------------------------------------------------------------
#---------------------------------------------------------------------


# Get nice plot background going
# (jet color for depth, land masked with grey)


# ---------------------------------------------
# Plot parameters
box_plot_mod = 2
box_plot_mod_island = 1
#box_plot_mod = 3

font_size_num_continent = 15
font_size_num_island = 12
font_size_labels = 17


# ---------------------------------------------


# I really don't know what "vmin" is doing here, but it seems to work (copied from a stack overflow post)

fig, ax = plt.subplots()

#mesh1 = ax.pcolormesh(lon_field,lat_field,mask_field,shading="nearest")
mesh1 = ax.pcolormesh(lon_field,lat_field,mask,shading="nearest")
ax.axis('image')

particle_numbers = particle_numbers_full

selected_files = [0] #W,S,Su,F
#selected_files = [0,2,4,7] #W,S,Su,F

#---------------------------------------------------------------------
# Testing, when working indent and make loop
#---------------------------------------------------------------------
#for tracking_output_file_pre in tracking_output_files:


dset = netCDF4.Dataset(tracking_output_file, 'r')

lon_pre = np.array(dset.variables['lon'])
lat_pre = np.array(dset.variables['lat'])
z_pre = np.array(np.abs(dset.variables['z']))
status = np.array(dset.variables['status'])
dset.close()


trajectory_mask = status == 0

for jj in range(np.shape(lon_pre)[0]):

    lon = lon_pre[jj,trajectory_mask[jj]]
    lat = lat_pre[jj,trajectory_mask[jj]]
    z = z_pre[jj,trajectory_mask[jj]]

    # plot initial locations - I can't yet figure out why some get deactivated.  Too close to the coast, stopped by the mask?
    ax.scatter(lon_pre[:,0], lat_pre[:,0],c='b',edgecolor='none')

    for kk in range(len(lon)-1):
        lon_seg = np.linspace(lon[kk],lon[kk+1],1000)
        lat_seg = np.linspace(lat[kk],lat[kk+1],1000)
        z_seg = np.linspace(z[kk],z[kk+1],1000)
    
        #z_segNorm = (z_seg - z_min)/z_range

        ax.scatter(lon_seg, lat_seg,c='g',edgecolor='none')
        ax.plot(lon[0],lat[0],'co')
        ax.plot(lon[-1],lat[-1],'ro')


#x_vecs_full = list(onshore_swim_component_x_map.flatten())
#y_vecs_full = list(onshore_swim_component_y_map.flatten())

#x_vecs_list = [x_vecs_full[ii] for ii in list(range(len(x_vecs_full))) if ii%10==0]
#y_vecs_list = [y_vecs_full[ii] for ii in list(range(len(y_vecs_full))) if ii%10==0]

#x_vecs_pre = np.array(x_vecs_list)
#y_vecs_pre = np.array(y_vecs_list)
#
#x_vecs = x_vecs_pre.reshape([len(x_vecs_pre),len(y_vecs_pre)])
#y_vecs = y_vecs_pre.reshape([len(x_vecs_pre),len(y_vecs_pre)])


x_vecs = np.zeros(np.shape(onshore_swim_component_x_map))
y_vecs = np.zeros(np.shape(onshore_swim_component_y_map))

for ii in range(np.shape(onshore_swim_component_x_map)[0]):
    for jj in range(np.shape(onshore_swim_component_x_map)[1]):
        if ii%10 == 0 and jj%10 == 0:
            x_vecs[ii,jj] = onshore_swim_component_x_map[ii,jj]
            y_vecs[ii,jj] = onshore_swim_component_y_map[ii,jj]

plt.quiver(lon_field,lat_field,x_vecs,y_vecs,color='r',scale=80)
#plt.quiver(lon_field,lat_field,onshore_swim_component_x_map,onshore_swim_component_y_map,color='r',scale=80)


plt.show()














