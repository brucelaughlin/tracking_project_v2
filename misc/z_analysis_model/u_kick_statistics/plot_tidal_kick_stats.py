


import pdb
import time

import datetime
import netCDF4
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as plt_path
from scipy.interpolate import interp1d
from geopy.distance import great_circle
import scipy.interpolate as spint
from os import listdir
from os.path import isfile, join
import sys
import argparse
from pathlib import Path
import os
import math

script_time_start = time.time()

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# PARAMETERS and CONSTANTS

# Opendrift swaps the particle and time dimensions relative to the tracking code that Patrick used


parser = argparse.ArgumentParser()
parser.add_argument("trackingdir", type=str)
#parser.add_argument("polygondex", type=int)
args = parser.parse_args()

tracking_output_dir = args.trackingdir
#study_polygon_index = args.polygondex

base_path = '/home/blaughli/tracking_project_v2/'
grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
grid_path_in = base_path + grid_directory + grid_file_in
dset = netCDF4.Dataset(grid_path_in, 'r')

lon_field = np.array(dset['lon_rho'])
lat_field = np.array(dset['lat_rho'])

dset.close

tracking_output_files = [f for f in listdir(tracking_output_dir) if isfile(join(tracking_output_dir,f))]
tracking_output_files.sort()

num_files = len(tracking_output_files)

# I always do this - is it bad practice?
tracking_output_dir = tracking_output_dir + "/"

#---------------------------------------------------------------------


#---------------------------------------------------------------------
# Get a lat/lon for a study... maybe also want to look at all kicks together?
#---------------------------------------------------------------------
tracking_output_file = tracking_output_dir + tracking_output_files[0]
dset = netCDF4.Dataset(tracking_output_file, 'r')
lons = dset.variables['lon'][:]
lats = dset.variables['lat'][:]
dset.close()

#study_polygon_index = 200
#num_releases_per_cell_per_seed = 5

#study_index = study_polygon_index * num_releases_per_cell_per_seed

#lon_study = lons[study_index,:][lons[study_index,:].mask == False][0]
#lat_study = lats[study_index,:][lats[study_index,:].mask == False][0]
#---------------------------------------------------------------------

file_number = 0

x_kicks = []
y_kicks = []

kick_speeds_all = []
bottom_depths_all = []

for tracking_output_file_pre in tracking_output_files:
    
    start_file_time = time.time()
   
    ###TESTING
#        if file_number > 0:
#            break

    print(tracking_output_file_pre)

    #tracking_output_file = os.path.abspath(tracking_output_file_pre)
    tracking_output_file = tracking_output_dir + tracking_output_file_pre

    dset = netCDF4.Dataset(tracking_output_file, 'r')

    #lons = np.array(dset.variables['lon'][:])
    #lats = np.array(dset.variables['lat'][:])
    #zs = np.array(dset.variables['z'][:])
    #kick_speeds = np.array(dset.variables['kick_speed'][:])
    #kick_angles = np.array(dset.variables['kick_angle'][:])
    
    lons = dset.variables['lon'][:]
    lats = dset.variables['lat'][:]
    zs = dset.variables['z'][:]
    kick_speeds = dset.variables['kick_speed'][:]
    kick_angles = dset.variables['kick_angle'][:]
    bottom_depths = dset.variables['sea_floor_depth_below_sea_level'][:]
    
    #lats = lats[lats.mask == False]
    
    #ocean_times = np.array(dset.variables['ocean_time'][:])
    ocean_times = np.array(dset.variables['time'][:])
    
    dset.close()


    num_bad = 0

    for particle_number in range(np.shape(lons)[0]):
        # Some particles have no kick - investigate
        if np.sum(np.shape(kick_speeds[particle_number,:][kick_speeds[particle_number,:].mask == False].data)) > 1:
            kick_speeds_all.append(kick_speeds[particle_number,:][kick_speeds[particle_number,:].mask == False][1])
            bottom_depths_all.append(bottom_depths[particle_number,:][bottom_depths[particle_number,:].mask == False][1])
        else:
            num_bad += 1
    


#        if lats[particle_number,:][lats[particle_number,:].mask == False][0] == lat_study and lons[particle_number,:][lons[particle_number,:].mask == False][0] == lon_study:
#       
#            # Guess there's no kick at first timestep?  That's the seeding timestep, ie position at time zero, right?
#            x_kicks.append(np.cos(kick_angles[particle_number,:][kick_angles[particle_number,:].mask == False][1]) * kick_speeds[particle_number,:][kick_speeds[particle_number,:].mask == False][1])
#            y_kicks.append(np.sin(kick_angles[particle_number,:][kick_angles[particle_number,:].mask == False][1]) * kick_speeds[particle_number,:][kick_speeds[particle_number,:].mask == False][1])


    break 

plt.scatter(kick_speeds_all,bottom_depths_all)
#plt.hist2d(x_kicks,y_kicks)
plt.colorbar()
plt.xlabel('x kick magnitude')
plt.ylabel('y kick magnitude')

plt.show()






