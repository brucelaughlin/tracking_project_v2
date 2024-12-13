# v1: trying to figure out why particles are being deactivated, specifically in the gradient swimming runs.  So, plotting
# over the swim mask, or can plot over rho mask

# https://stackoverflow.com/questions/20165169/change-colour-of-curve-according-to-its-y-value-in-matplotlib

# Input Files
#---------------------------------------------------------------------

#tracking_output_dir_pre = "drift_150_swim_test_v2_zeroVel_grad/"
#swim_data_file = '/home/blaughli/tracking_project/practice/onshore_swim_work/z_production/swim_data.p'

tracking_output_dir_pre = "drift_150_swim_test_v2_zeroVel/"
swim_data_file = '/home/blaughli/tracking_project/practice/onshore_swim_work/z_production/swim_data_simple.p'

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
tracking_output_base = "/data01/blaughli/tracking_project_output/"
#tracking_output_base = "/data03/blaughli/tracking_project_output/"
tracking_output_dir = tracking_output_base + tracking_output_dir_pre

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


tracking_output_files = [f for f in listdir(tracking_output_dir) if isfile(join(tracking_output_dir,f))]
tracking_output_files.sort()


file = open(swim_data_file,'rb')
mask,mask_flat,coord_array,onshore_swim_component_x,onshore_swim_component_y,onshore_swim_component_x_map,onshore_swim_component_y_map = pickle.load(file)
#mask,mask_flat,coord_array,onshore_swim_component_x,onshore_swim_component_y,neg_grad_norm_map_x,neg_grad_norm_map_y = pickle.load(file)
#mask,mask_flat,coord_array,onshore_swim_component_x,onshore_swim_component_y = pickle.load(file)
#mask_flat,coord_array,onshore_swim_component_x,onshore_swim_component_y = pickle.load(file)
file.close()


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

tFiles = [tracking_output_files[jj] for jj in selected_files]

# First, determine absolue min/max and range of Z, for consistent plotting and colorbar
#z_range = 0 
#z_min = 999999
#z_max = 0
#for ii in range(len(tFiles)):
#    tracking_output_file_pre = tFiles[ii]
#    tracking_output_file = tracking_output_dir + tracking_output_file_pre
#    dset = netCDF4.Dataset(tracking_output_file, 'r')
#    z_pre = np.abs(dset.variables['z'][particle_numbers])
#    dset.close()
#    z_max_pre = np.max(z_pre)
#    z_min_pre = np.min(z_pre)
#    if z_min_pre < z_min:
#        z_min = z_min_pre
#    if z_max_pre > z_max:
#        z_max = z_max_pre
#z_range = z_max - z_min


for ii in range(len(tFiles)):

    tracking_output_file_pre = tFiles[ii]
    tracking_output_file = tracking_output_dir + tracking_output_file_pre

    dset = netCDF4.Dataset(tracking_output_file, 'r')

    lon_pre = np.array(dset.variables['lon'])
    lat_pre = np.array(dset.variables['lat'])
    z_pre = np.array(np.abs(dset.variables['z']))
    status = np.array(dset.variables['status'])
    #lon_pre = dset.variables['lon'][particle_numbers]
    #lat_pre = dset.variables['lat'][particle_numbers]
    #z_pre = np.abs(dset.variables['z'][particle_numbers])
    #status = dset.variables['status'][particle_numbers]
    dset.close()

    #color_steps = np.linspace(0,1,1000)

    trajectory_mask = status == 0

    for jj in range(np.shape(lon_pre)[0]):
    #for jj in range(len(particle_numbers)):
#---------------------------------------------------------------------
#        jj = 0

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
            #ax.scatter(lon_seg, lat_seg, c = cm.jet(z_segNorm), edgecolor='none')
            ax.plot(lon[0],lat[0],'co')
            ax.plot(lon[-1],lat[-1],'ro')

plt.quiver(lon_field,lat_field,onshore_swim_component_x_map.T,onshore_swim_component_y_map.T,color='r',scale=80)

#N = 20
#cmap = plt.get_cmap('jet',N)
#
#cbar_left = 0.8
#cbar_bot = 0.2
#cbar_width = 0.02
#cbar_height = 0.5
#
##sub_ax = plt.axes([0.96, 0.55, 0.02, 0.3]) 
#sub_ax = plt.axes([cbar_left,cbar_bot,cbar_width,cbar_height]) 
#
##norm = mpl.colors.Normalize(vmin=z_min,vmax=z_max)
#norm = plt.Normalize(vmin=z_min,vmax=z_max)
#sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#sm.set_array([])
##cbar = plt.colorbar(sm, ticks=np.linspace(z_min,z_max,N+1), cax = sub_ax, boundaries=np.arange(-0.05,2.1,.1))
#cbar = plt.colorbar(sm, ticks=np.linspace(z_min,z_max,N+1), cax = sub_ax)
##cbar = plt.colorbar(sm, ticks=np.linspace(z_min,z_max,N), cax = sub_ax)
#
#cbar_label = "Depth (m)"
#cbar_fontSize = 10
#
#cbar.ax.set_ylabel(cbar_label, fontsize = cbar_fontSize)



plt.show()














