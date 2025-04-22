
# oldConCalc: Trying to use the old "trajectory by trajectory" way of calculating connectivity - seems to be much faster

# cleanButSlower:  I tried to "vectorize" this calculation more cleanly, but my first efforts (see "_dumb" in x_old)
# were misguided because I was having to process every single timestep rather than just those in a given PLD (this made the calculation take forever).
# Now, I've taken what I liked from those and used ideas from the previous working versions, and the code is cleaner but takes ~twice as long
# as the original kludgy version.  So, something I did originally was relatively fast.  

# In the next version, I'd like to use what works from both versions to make a clean production version...


script_version = "cleanButSlower_oldConCalc"

# Note that "status" is 0 when the particle is active, and a large magnitude negative
# number when not.  (strange!)

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
parser.add_argument("pdrakeswitch", nargs='?', type=str)
args = parser.parse_args()

tracking_output_dir = args.trackingdir
pdrake_switch = args.pdrakeswitch


save_output_directory = '/home/blaughli/tracking_project_v2/processing/settlement/z_output'
polygon_file_path = '/home/blaughli/tracking_project_v2/input_files/wc15.0_06.0km_036km2_settling_polygons.txt'



# For my files, the particle dimension is the first dimension (ie 0):
particle_dimension = 0
time_dimension = 1
# For Patrick's files, the particle dimension is the second dimension (ie 1):
if pdrake_switch != None:
    particle_dimension = 1
    time_dimension = 0


#---------------------------------------------------------------------
# Define PLDs
#---------------------------------------------------------------------
# (Using PLDs for the species listed in Mallarie Yeager's report)


# MY NEW ALGORITHM REQUIRES THE ADDITION OF THE "[0,1]" PLD TO THE BEGINNING OF THE ARRAY

# (THIS FIRST PLD IS USED TO DETERMINE INITIAL LOCATION"
pld_seed = [0,1]
## KELP BASS
pld_kelp_bass = [20,29]
## CALIFORNIA SHEEPHEAD
pld_ca_sheephead = [30,59]
## KELP ROCKFISH
pld_kelp_rockfish = [60,89]
## BLUE ROCKFISH, BLACK ROCKFISH
pld_blue_black_rockfish = [90, 149]

#pld_test = [170,174]
#pld_test = [45,49]
pld_test = [8,9]
#pld_test = [5,9]

pld_list_array = np.array([[0,1],[45,49]])
#pld_list_array = np.array([[45,49]])
#pld_list_array = np.array([[5,6],[10,11],[15,17],[20,22],[30,33],[45,49],[60,65],[90,98],[120,131]])


num_pld_with_0 = len(pld_list_array)
num_pld = num_pld_with_0 - 1
#num_pld = len(pld_list_array)
#num_pld_with_0 = num_pld + 1



#---------------------------------------------------------------------
#---------------------------------------------------------------------

# This should point to my own polygons!!  Still need to put those into a CSV txt file using patrick's file as a template!


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


#---------------------------------------------------------------------
# Use a text csv file of polygon vertices to create the boxes
# These files must have a consistent format for their lines:
#   cell #, vertex #, lat, lon
#   001, 01, 30.050000, -115.816661
# Cell numbers in the file must always start at 1!!!!!

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
#---------------------------------------------------------------------

print('num polygons: {}'.format(num_polygons))


extrema_lon = [np.min(lon_field),np.max(lon_field)]
extrema_lat = [np.min(lat_field),np.max(lat_field)]


#---------------------------------------------------------------------

# number of pdfs per feild (1 annual, 4 seasonal = 5 total)
# JUST MAKE THE ANNUAL FILE
n_pdfs = 1
#n_pdfs = 5


# Need to have a way to see if I'm doing this right... without double counting errors, etc
tracking_output_file = os.path.join(tracking_output_dir,tracking_output_files[0])
#tracking_output_file = tracking_output_dir + tracking_output_files[0]
dset = netCDF4.Dataset(tracking_output_file, 'r')

if time_dimension == 0:
    ocean_time = dset.variables['ocean_time'][:]
else:
    ocean_time = dset.variables['time'][:]
lon = dset.variables['lon'][:]
dset.close()
num_trajectories = np.shape(lon)[particle_dimension]
num_files = len(tracking_output_files)

# I don't think <num_trajectories> is necessarily fixed (final seeding in run might have fewer seeds than other files)
counter_array=np.zeros((num_trajectories,num_files))


print('NUMBER OF PARTICLES: {}'.format(num_trajectories))


timesteps_per_day = int(1/((ocean_time[1]-ocean_time[0])/86400))
#timesteps_per_day = int(1/((ocean_ocean_time[1]-ocean_ocean_time[0])/86400))

# For now, exit if output is less than daily...
if timesteps_per_day < 1:
    raise RuntimeError('This code is dumb and currently breaks if the save timestep is > 1 day.  Should be an easy fix')


pld_list_array_timesteps = pld_list_array/timesteps_per_day
pld_list_array_timesteps = pld_list_array_timesteps.astype(int)

#----------------------------------------------------------
# Prepare exposure information
#----------------------------------------------------------

# For the histogram of average temperature experienced - just estimating range, since I don't know it without processing
#---------------------
# Make sure these match (0.1 = 1, 0.01 = 2, etc)
T_step = 0.1
n_decimals_round = 1
T_scale_factor = int(1/T_step)
#---------------------
T_min = 0
T_max = 30
n_T_bins = len(np.arange(T_min,T_max+1,T_step))

# Looking at Jerome's files:
#    float O2(ocean_time,s_rho,eta_rho,xi_rho) ;
#      O2:long_name = "time-averaged dissolved O2 concentration" ;
#      O2:units = "millimole_O2 meter-3" ;

# desired units: mg/L
molarMassO2 = 31.999 # g/mol
conversion_factor = molarMassO2/1000  #worked this out on paper

# Limit below which we care about exposure for O2
O2_limit_list = [2.2,3.1,4.1,6]
pH_limit_list = [7.5,7.6,7.7,7.8,7.9,8,8.1,8.2,8.3]



#---------------------------------------------------------------------
# START THE MAIN LOOP!!!
#---------------------------------------------------------------------

        
#---------------------------------------------------------------------
# Create lists to store statistics for full run and seasonal subsets
#---------------------------------------------------------------------

# store the number of particles released from each box
release_counts_per_cell = np.zeros((n_pdfs,num_polygons), dtype=int)




# Connectivity (release box number vs settlement box number)
pdf_arrays_connectivity = np.zeros((num_pld,n_pdfs,num_polygons,num_polygons))
#pdf_arrays_connectivity = np.zeros((num_pld,n_pdfs,num_polygons,num_polygons))
#pdf_arrays_connectivity = np.zeros((n_pdfs,num_polygons,num_polygons))



# Time after PLD until settlement (saving only release location) (release box number vs settlement time)
#pdf_arrays_settleTime = np.zeros((n_pdfs,num_polygons,timesteps_settlement_window))

# Number of days eposed to DO levels below 2.2 (saving only settlement location) (release box number vs exposure time)
# note: Think the time dimension needs to be one bigger than the number of possible timesteps (ie need to include an option for "zero") 
#pdf_arrays_O2 = np.zeros((len(O2_limit_list),n_pdfs,num_polygons,timesteps_full_run))

# Same idea for pH
#pdf_arrays_pH = np.zeros((len(pH_limit_list),n_pdfs,num_polygons,timesteps_full_run))

# Histogram of average temperature experienced - just estimating range, since I don't know it without processing
#pdf_arrays_T = np.zeros((n_pdfs,num_polygons,n_T_bins))
   

# Must track the "fake" local settlers, ie particles that never move!!!
#stationary_settlers_array = np.zeros((num_polygons,num_files))

#empty_trajectory_count_array = np.zeros(num_files)


#---------------------------------------------------------------------
#---------------------------------------------------------------------



file_number = 0

for tracking_output_file_pre in tracking_output_files:

    file_time_start = time.time()

    ###TESTING
    #if file_number > 1:
    #if file_number > 0:
        #break

    print(tracking_output_file_pre)

        
    tracking_output_file = tracking_output_dir + tracking_output_file_pre

    dset = netCDF4.Dataset(tracking_output_file, 'r')

    if time_dimension == 0:
        ocean_time = dset.variables['ocean_time'][:]
    else:
        ocean_time = dset.variables['time'][:]

    lon_all = dset.variables['lon'][:]
    lat_all = dset.variables['lat'][:]
    
    dset.close()


    # Store the total number of particles                                                                                                                     
    num_trajectories = np.shape(lon_all)[particle_dimension]                                                                                                     
    num_timesteps = np.shape(lon_all)[time_dimension]                                                                                                         
    
    lons = np.full([num_trajectories,num_timesteps], np.nan)
    lats = np.full([num_trajectories,num_timesteps], np.nan)
       
    for trajectory_dex in range(num_trajectories):                                                                                                                 
                                                                                                                                                              
        if time_dimension == 0:                                                                                                                               
            trajectory_lon = lon_all[:,trajectory_dex]
            trajectory_lat = lat_all[:,trajectory_dex]
        else:                                                                                                                                                 
            trajectory_lon = lon_all[trajectory_dex,:]
            trajectory_lat = lat_all[trajectory_dex,:]
        trajectory_lon = trajectory_lon[np.logical_not(trajectory_lon.mask)].data                                                                               
        trajectory_lat = trajectory_lat[np.logical_not(trajectory_lat.mask)].data                                                                               
        lons[trajectory_dex,0:len(trajectory_lon)] = trajectory_lon
        lats[trajectory_dex,0:len(trajectory_lat)] = trajectory_lat
                                                                                                                                                              


    # For keeping track of settlement data
#    pld_indices = np.zeros(num_trajectories, dtype=int)
#    pld_masks = np.ones((num_trajectories,num_pld_with_0), dtype=bool)
    pld_data = np.zeros((num_trajectories,num_pld_with_0))
#    still_drifting_indices = np.ones(num_trajectories, dtype=bool)


    for pld_dex in range(np.shape(pld_list_array)[0]):
        pld_settle_mask = np.ones(num_trajectories, dtype=bool)
        for tt in range(pld_list_array_timesteps[pld_dex,0],pld_list_array_timesteps[pld_dex,1]): 

            print("file {}/{}, timestep {}/{}".format(file_number+1, num_files, tt, pld_list_array_timesteps[pld_dex,1] - 1))

            if time_dimension == 0:
                lon_tt = lons[tt,:]
                lat_tt = lats[tt,:]
            else:
                lon_tt = lons[:,tt]
                lat_tt = lats[:,tt]

                # Exposure variables 
                #O2_all = dset.variables['oxygen'][:,tt]
                #pH_all = dset.variables['pH'][:,tt]
                #temp_all = dset.variables['sea_water_temperature'][:,tt]
                #O2_all *= conversion_factor
            # add settlement time storage
            #settlement_times = np.full((num_trajectories), np.nan)

            points_lon_lat = np.array([lon_tt,lat_tt]).T # originally built nx2 matrices, so i guess i'll stick with two-column matrices
            timestep_timer_start = time.time()
            for polygon_dex in range(num_polygons):
                path = plt_path.Path(polygon_vertex_list[polygon_dex])
                internal_particles_indices = path.contains_points(points_lon_lat)
           
                settler_indices = np.logical_and(internal_particles_indices, pld_settle_mask)
        
                if np.sum(settler_indices) == 0:
                    print('No settlers at current timestep')
                    continue
                
                pld_data[settler_indices, pld_dex] = polygon_dex + 1
                pld_settle_mask[settler_indices] = False
                #pld_settle_mask[internal_particles_indices] = False

            timestep_timer_end = time.time()
            print(f'Timestep loop took {(timestep_timer_end-timestep_timer_start)/60:.03f} minutes')
        


## THESE SHOULD BE MULTIPLE FUNCTION CALLS HANDLING THE INDIVIDUAL TASKS (CALCULATE SETTLEMENT LOCATIONS, CREATE CONNECTIVITY MATRIX, ETC)

        #if any(pld_data[:,0] < 0):
        #if pld_data[] < 0:
        #if release_boxes[trajectory_dex] < 0:
        #raise RuntimeError('{} particles had negative release cell indices}'.format(sum(erronious_seeding_indices)))
            #raise RuntimeError('particle {} has release cell dex: {}'.format(trajectory_dex,release_boxes[trajectory_dex]))
        # --------------------------------------------------
        # --------------------------------------------------
        
      
   
    # NOW FOR CONNECTIVITY MATRIX (MAKE THIS A FUNCTION)

    erronious_seeding_indices = pld_data[:,0] < 0
    box_zero_seeding_indices = pld_data[:,0] == 0

    print(f'{sum(erronious_seeding_indices)} particles had negative release cell indices')
    print(f"{sum(box_zero_seeding_indices)} particles were recorded released in 'box zero', the phantom box")
       


    # Patrick only has the annual array
    array_dex_list = [0]
  
    con_time_start = time.time()

    release_counts_per_cell = np.zeros(num_polygons)
    for polygon_dex in range(num_polygons):
        release_counts_per_cell[polygon_dex] += np.sum(pld_data[:,0] == polygon_dex + 1) 

   
    # pld_indices = np.zeros(num_trajectories)  
    # pld_masks = np.ones(num_trajectories,num_pld+1, dtype=bool)
    # pld_data = np.zeros(num_trajectories,num_pld+1)


    for pld_dex in range(1,num_pld_with_0):
        pld_dex_discounting_release_pld = pld_dex - 1
        for array_dex in array_dex_list:
            for trajectory_dex in range(np.shape(lons)[0]):
                #if np.logical_not(math.isnan(settlement_boxes[trajectory_dex])):
                    
                pdf_arrays_connectivity[pld_dex_discounting_release_pld,array_dex,int(pld_data[trajectory_dex,0])-1,int(pld_data[trajectory_dex,pld_dex])-1] += 1

                #pdf_arrays_connectivity[array_dex,int(release_boxes[trajectory_dex])-1,int(settlement_boxes[trajectory_dex])-1] += 1 #from older script







    #for pld_dex in range(1,num_pld_with_0):
    #    pld_dex_discounting_release_pld = pld_dex - 1
    #    for array_dex in array_dex_list:
    #        for ii in range(num_polygons):
    #            for jj in range(num_polygons):
    #                pdf_arrays_connectivity[pld_dex_discounting_release_pld,array_dex,ii,jj] += np.sum(np.logical_and( pld_data[:,0] == ii + 1, pld_data[:,pld_dex] == jj + 1 ))   
                
                
            
    con_time_end = time.time()
    
    print(f'Connectivity loop took {(con_time_end-con_time_start)/60:.03f} minutes')

    file_number += 1

    file_time_end = time.time()
    print(f'File loop took {(file_time_end-file_time_start)/60:.03f} minutes')

print(f'Total processing time: {(file_time_end-script_time_start)/60:.03f} minutes')





#print('Number of stationary settlers: {}'.format(num_stationary_settlers))

d = {}
d['release_counts_per_cell'] = release_counts_per_cell
d['pdf_arrays_connectivity'] = pdf_arrays_connectivity
#d['pdf_arrays_settleTime'] = pdf_arrays_settleTime
#d['counter_array'] = counter_array
#d['stationary_settlers_array'] = stationary_settlers_array
#d['empty_trajectory_count_array'] = empty_trajectory_count_array
d['pld_days'] = pld_list_array[pld_dex]

#if time_dimension == 1:
#    d['pdf_arrays_T'] = pdf_arrays_T
#    d['pdf_arrays_O2'] = pdf_arrays_O2
#    d['pdf_arrays_pH'] = pdf_arrays_pH
#    d['O2_limit_list'] = O2_limit_list
#    d['pH_limit_list'] = pH_limit_list


save_output_file_name_pre = "connectivity_data_"
#save_output_file_name_pre = "binned_data_seasonal_allReleases_baseYear_{}_".format(base_year)
save_output_file_name = save_output_file_name_pre + tracking_output_dir.split('/')[-2]
save_output_full_path = os.path.join(save_output_directory, save_output_file_name + "_pld_{}_{}_version_{}".format(pld_list_array[pld_dex,0],pld_list_array[pld_dex,1], script_version))

np.savez(save_output_full_path, **d)







