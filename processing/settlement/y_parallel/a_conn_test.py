# Copied form "basicVars" version.  Now parallelizing


#script_version = "basicVars"
script_version = "singleFileParallel"


#save_output_directory = '/home/blaughli/tracking_project_v2/processing/settlement/z_output'

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

#parser.add_argument("--configfile", type=str, help='help yourself')

parser = argparse.ArgumentParser()
parser.add_argument("--trackingfile", type=str)
#parser.add_argument("trackingdir", type=str)
parser.add_argument("--pdrakefileswitch", type=int)
parser.add_argument("--pdrakepolygonswitch", type=int)
#parser.add_argument("pdrakeswitch", nargs='?', type=str)
args = parser.parse_args()

tracking_file = args.trackingfile
#tracking_dir = args.trackingfile
pdrake_file_switch = args.pdrakefileswitch
pdrake_polygon_switch = args.pdrakepolygonswitch

tracking_dir = os.path.dirname(tracking_file)


# For my files, the particle dimension is the first dimension (ie 0):
particle_dimension = 0
time_dimension = 1
# For Patrick's files, the particle dimension is the second dimension (ie 1):
if pdrake_file_switch != 0:
    particle_dimension = 1
    time_dimension = 0

if pdrake_polygon_switch != 0:
    # pDrake Cells
    polygon_file_path = '/home/blaughli/tracking_project_v2/input_files/wc15.0_06.0km_036km2_settling_polygons.txt'
    polygon_string = "pDrake_polygons"
else:
    # bLaughli Cells
    polygon_file_path = '/home/blaughli/tracking_project_v2/misc/wc15n_300km2_settling_polygons.txt'
    polygon_string = "bLaughli_polygons"




#---------------------------------------------------------------------
# Define PLDs
#---------------------------------------------------------------------
# (Using PLDs for the species listed in Mallarie Yeager's report)


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

pld_array = np.array([[5,6],[10,11],[15,17],[20,22],[30,33],[45,49],[60,65],[90,98],[120,131]])
#pld_array=np.array([pld_kelp_bass,pld_ca_sheephead,pld_kelp_rockfish,pld_blue_black_rockfish,pld_test])
#pld_array=np.array([pld_kelp_bass,pld_ca_sheephead,pld_kelp_rockfish,pld_blue_black_rockfish])

# choose one pld to use, while still testing
#pld_chosen_dex = len(pld_array)-1
#pld_chosen_dex = 2
pld_chosen_dex = 5

#---------------------------------------------------------------------
#---------------------------------------------------------------------

base_path = '/home/blaughli/tracking_project_v2/'
grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
grid_path_in = base_path + grid_directory + grid_file_in
dset = netCDF4.Dataset(grid_path_in, 'r')

lon_field = np.array(dset['lon_rho'])
lat_field = np.array(dset['lat_rho'])

dset.close

bounding_boxes_base = base_path + 'practice/bounding_boxes/create_boxes/'
bounding_boxes_continent_dir = bounding_boxes_base + 'continent/z_output/'
bounding_boxes_islands_dir = bounding_boxes_base + 'modify_islands/z_output/'


# I always do this - is it bad practice?
tracking_dir = tracking_dir + "/"

#save_output_directory = '/home/blaughli/tracking_project_v2/processing/settlement/z_output'

save_output_directory = tracking_dir + 'v_connHist_files'

#if not os.path.exists(save_output_directory):
#    os.mkdirs(save_output_directory)
    
Path(save_output_directory).mkdir(parents=True, exist_ok=True)

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

# number of pdfs per feild (4 seasonal and 1 overall = 5)
# JUST MAKE THE ANNUAL FILE
n_pdfs = 1
#n_pdfs = 5


dset = netCDF4.Dataset(tracking_file, 'r')

if time_dimension == 0:
    ocean_time = dset.variables['ocean_time'][:]
else:
    ocean_time = dset.variables['time'][:]
lon = dset.variables['lon'][:]
dset.close()
num_particles = np.shape(lon)[particle_dimension]
counter_array=np.zeros((num_particles))


print('NUMBER OF PARTICLES: {}'.format(num_particles))


timesteps_per_day = int(1/((ocean_time[1]-ocean_time[0])/86400))
#timesteps_per_day = int(1/((ocean_ocean_time[1]-ocean_ocean_time[0])/86400))

# For now, exit if output is less than daily...
if timesteps_per_day < 1:
    raise RuntimeError('This code is dumb and currently breaks if the save timestep is > 1 day.  Should be an easy fix')


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


# Delete this, just index by length
dummy_value = 9999

#---------------------------------------------------------------------
# START THE MAIN LOOP!!!
#---------------------------------------------------------------------

for pld_dex in range(len(pld_array)):
        
    ###TESTING
    # test using the chosen pld
    if pld_dex != pld_chosen_dex:
        continue
        #break

    first_settlement_day = pld_array[pld_dex,0]
    last_settlement_day = pld_array[pld_dex,1]

    print('pld: {}-{}'.format(first_settlement_day,last_settlement_day))

    # THIS is the v5 adjustment: need to only use data from the pld, which begins after "first_settlement_day" and ends after "last_settlement_day"
    pld_length_days = last_settlement_day - first_settlement_day + 1


    timesteps_settlement_window = pld_length_days * timesteps_per_day
    timesteps_full_run = (last_settlement_day+1) * timesteps_per_day + 1 # why isn't this just taken from the dimensions of the files?


    #first_settle_dex = first_settlement_day * timesteps_per_day + 1
    first_settle_dex = first_settlement_day * timesteps_per_day  # Why was I adding 1 here????
    last_settle_dex = (last_settlement_day+1) * timesteps_per_day # And here????
    #last_settle_dex = (last_settlement_day+1) * timesteps_per_day + 1


    #---------------------------------------------------------------------
    # Create lists to store statistics for full run and seasonal subsets
    #---------------------------------------------------------------------

    # store the number of particles released from each box
    release_counts_per_cell = np.zeros((n_pdfs,num_polygons), dtype=int)


    # Connectivity (release box number vs settlement box number)
    pdf_arrays_connectivity = np.zeros((n_pdfs,num_polygons,num_polygons))

    # Time after PLD until settlement (saving only release location) (release box number vs settlement time)
    pdf_arrays_settleTime = np.zeros((n_pdfs,num_polygons,timesteps_settlement_window))

    # Number of days eposed to DO levels below 2.2 (saving only settlement location) (release box number vs exposure time)
    # note: Think the time dimension needs to be one bigger than the number of possible timesteps (ie need to include an option for "zero") 
    pdf_arrays_O2 = np.zeros((len(O2_limit_list),n_pdfs,num_polygons,timesteps_full_run))

    # Same idea for pH
    pdf_arrays_pH = np.zeros((len(pH_limit_list),n_pdfs,num_polygons,timesteps_full_run))

    # Histogram of average temperature experienced - just estimating range, since I don't know it without processing
    pdf_arrays_T = np.zeros((n_pdfs,num_polygons,n_T_bins))
       
    
    # Must track the "fake" local settlers, ie particles that never move!!!
    stationary_settlers_array = np.zeros((num_polygons,1))

    empty_trajectory_count_array = np.zeros(1)


    #---------------------------------------------------------------------
    #---------------------------------------------------------------------

    file_number = 0

    start_file_time = time.time()
   
    ###TESTING
#        if file_number > 0:
#            break

    print(tracking_file)

    #tracking_file_full_path = tracking_output_dir + tracking_file

    dset = netCDF4.Dataset(tracking_file, 'r')
    #dset = netCDF4.Dataset(tracking_file_full_path, 'r')

    lon_all = dset.variables['lon'][:]
    lat_all = dset.variables['lat'][:]
    #z_all = dset.variables['z'][:]
    if time_dimension == 0:
        ocean_time = dset.variables['ocean_time'][:]
    else:
        ocean_time = dset.variables['time'][:]
        # Exposure variables 
        #O2_all = dset.variables['oxygen'][:]
        #pH_all = dset.variables['pH'][:]
        temp_all = dset.variables['sea_water_temperature'][:]
        
        #O2_all *= conversion_factor

    dset.close()








    stationary_settlers_per_box = np.zeros((num_polygons,1))

    # Store the total number of particles
    num_particles = np.shape(lon_all)[particle_dimension]
    num_timesteps = np.shape(lon_all)[time_dimension]
     
    # Why not just make one array of all of the data, and then reference it by time index???
    if time_dimension == 0:
        particles_lon_all=np.full((timesteps_full_run,num_particles), np.nan)
        particles_lat_all=np.full((timesteps_full_run,num_particles), np.nan)
        #particles_O2_all=np.zeros([timesteps_full_run,num_particles])
        #particles_pH_all=np.zeros([timesteps_full_run,num_particles])
        particles_T_all=np.zeros([timesteps_full_run,num_particles])
    else: 
        particles_lon_all=np.full((num_particles,timesteps_full_run), np.nan)
        particles_lat_all=np.full((num_particles,timesteps_full_run), np.nan)
        #particles_O2_all=np.zeros([num_particles,timesteps_full_run])
        #particles_pH_all=np.zeros([num_particles,timesteps_full_run])
        particles_T_all=np.zeros([num_particles,timesteps_full_run])
    
    for particle_num in range(num_particles):
        
        if time_dimension == 0:
            particle_lon = lon_all[:,particle_num]
            particle_lon = particle_lon[np.logical_not(particle_lon.mask)].data
            particle_lon = particle_lon[0:timesteps_full_run]
            particle_lat = lat_all[:,particle_num]
            particle_lat = particle_lat[np.logical_not(particle_lat.mask)].data
            particle_lat = particle_lat[0:timesteps_full_run]
        else:
            particle_lon = lon_all[particle_num,:]
            particle_lon = particle_lon[np.logical_not(particle_lon.mask)].data
            particle_lon = particle_lon[0:timesteps_full_run]
            particle_lat = lat_all[particle_num,:]
            particle_lat = particle_lat[np.logical_not(particle_lat.mask)].data
            particle_lat = particle_lat[0:timesteps_full_run]
        
            #particle_O2 = O2_all[particle_num,:]
            #particle_O2 = particle_O2[np.logical_not(particle_O2.mask)].data
            #particle_O2 = particle_O2[0:timesteps_full_run]
            #particle_pH = pH_all[particle_num,:]
            #particle_pH = particle_pH[np.logical_not(particle_pH.mask)].data
            #particle_pH = particle_pH[0:timesteps_full_run]
            particle_T = temp_all[particle_num,:]
            particle_T = particle_T[np.logical_not(particle_T.mask)].data
            particle_T = particle_T[0:timesteps_full_run]
        
        if time_dimension == 0:
            particles_lon_all[0:len(particle_lon),particle_num] = particle_lon
            particles_lat_all[0:len(particle_lat),particle_num] = particle_lat
        else:
            particles_lon_all[particle_num,0:len(particle_lon)] = particle_lon
            particles_lat_all[particle_num,0:len(particle_lat)] = particle_lat
            #particles_O2_all[particle_num,:] = np.pad(particle_O2,(0,timesteps_full_run-len(particle_O2)), 'constant',constant_values=(dummy_value))
            #particles_pH_all[particle_num,:] = np.pad(particle_pH,(0,timesteps_full_run-len(particle_pH)), 'constant',constant_values=(dummy_value))
            particles_T_all[particle_num,:] = np.pad(particle_T,(0,timesteps_full_run-len(particle_T)), 'constant',constant_values=(dummy_value))

    # -------------
    # New idea for exposure, given the difficulty of removing stationary settlers AFTER histograms have been calculated.
    # Why not just threshold breaches beforehand, and then ... keep track somehow of which points have been added which should be removed...
    
    #if time_dimension == 1:
    #    O2_threshold_breach_array = np.zeros((len(O2_limit_list),np.shape(particles_O2_all)[0],np.shape(particles_O2_all)[1]), dtype=bool)
    #    pH_threshold_breach_array = np.zeros((len(pH_limit_list),np.shape(particles_pH_all)[0],np.shape(particles_pH_all)[1]), dtype=bool)

    #pdb.set_trace()
        #for limit_dex in range(len(O2_limit_list)):
        #    O2_threshold_breach_array[limit_dex,:,:] = particles_O2_all < O2_limit_list[limit_dex]
        #for limit_dex in range(len(pH_limit_list)):
        #    pH_threshold_breach_array[limit_dex,:,:] = particles_pH_all < pH_limit_list[limit_dex]

        #particles_O2_all = np.where(particles_O2_all == dummy_value, np.nan, particles_O2_all)
        #particles_O2_all = np.ma.array(particles_O2_all, mask = np.isnan(particles_O2_all))
        
        #particles_pH_all = np.where(particles_pH_all == dummy_value, np.nan, particles_pH_all)
        #particles_pH_all = np.ma.array(particles_pH_all, mask = np.isnan(particles_pH_all))
        
        particles_T_all = np.where(particles_T_all == dummy_value, np.nan, particles_T_all)
        particles_T_all = np.ma.array(particles_T_all, mask = np.isnan(particles_T_all))
        
        #O2_threshold_breach_array = np.where(O2_threshold_breach_array == dummy_value, np.nan, O2_threshold_breach_array)
        #O2_threshold_breach_array = np.ma.array(O2_threshold_breach_array, mask = np.isnan(O2_threshold_breach_array))
        
        #pH_threshold_breach_array = np.where(pH_threshold_breach_array == dummy_value, np.nan, pH_threshold_breach_array)
        #pH_threshold_breach_array = np.ma.array(pH_threshold_breach_array, mask = np.isnan(pH_threshold_breach_array))
        
        # Exposure statistics
        #particle_arrays_O2 = np.zeros((num_particles,len(O2_limit_list)), dtype=int)
        #particle_arrays_pH = np.zeros((num_particles,len(pH_limit_list)), dtype=int)
        particle_array_T = np.zeros(num_particles)
        particle_array_driftTime = np.zeros(num_particles)


    # LOOK AT NP.MAWHERE - ACCOMPLISH IN ONE STEP
    # or np.maand - may be an "and" for masked arrays, which takes masks into account
    particles_lon_all = np.where(particles_lon_all == dummy_value, np.nan, particles_lon_all)
    particles_lon_all = np.ma.array(particles_lon_all, mask = np.isnan(particles_lon_all))
    
    particles_lat_all = np.where(particles_lat_all == dummy_value, np.nan, particles_lat_all)
    particles_lat_all = np.ma.array(particles_lat_all, mask = np.isnan(particles_lat_all))
    


    # Create array for tracking drift time, for use in computing average T along trajectory 
    particle_drift_times = np.zeros(num_particles)

    # prepare lists to hold the starting/ending box numbers for each particle
    release_boxes = np.zeros(num_particles, dtype=int)
    settlement_boxes = np.zeros(num_particles, dtype=int)
    # add settlement time storage
    settlement_times = np.zeros(num_particles, dtype=int)


    #---------------------------------------------------------------------
    # FIRST MUST DETERMINE STARTING LOCATIONS
    #---------------------------------------------------------------------
    
    # <points_lon_lat> is used here and in settlement; for starting locations, its value is fixed

    points_lon_lat = np.zeros((num_particles,2))
    if time_dimension == 0:
        points_lon_lat[:,0] = particles_lon_all[0,:]
        points_lon_lat[:,1] = particles_lat_all[0,:]
    else:
        points_lon_lat[:,0] = particles_lon_all[:,0]
        points_lon_lat[:,1] = particles_lat_all[:,0]

    particle_safety_mask = np.ones(num_particles, dtype=bool)
    
    for polygon_dex in range(num_polygons):
        path = plt_path.Path(polygon_vertex_list[polygon_dex])
        particles_inside_flags = path.contains_points(points_lon_lat)
        
        release_boxes[np.logical_and(particles_inside_flags,particle_safety_mask)] = polygon_dex + 1
        particle_safety_mask[particles_inside_flags] = False

    
    # Need a mask to prevent processing of particles with no release location (I have no idea how that happens, by the way..???!)
    empty_release_cells = np.setdiff1d(np.arange(num_polygons)+1,release_boxes)

    print("--------------------------------------------")
    print(f"{sum(empty_release_cells)} Empty release cells")
    print("--------------------------------------------")

    
    # NOW DETERMINE SETTLEMENT LOCATIONS!

    # create the "safety mask" that I'll use to make sure I only store the first box entered during the settlement window
    particle_safety_mask = np.ones(num_particles, dtype=bool)

    stationary_settler_array = np.ones((num_particles,num_polygons), dtype=bool) 

    for time_dex in range(first_settle_dex,last_settle_dex):

        # TESTING
        #if time_dex > first_settle_dex + 1:
        #if time_dex > first_settle_dex:
        #    break

        print("file {}/{}, timestep {}/{}".format(file_number+1, 1, time_dex,last_settle_dex-1))

        points_lon_lat = np.zeros((num_particles,2))
        if time_dimension == 0:
            points_lon_lat[:,0] = particles_lon_all[time_dex,:]
            points_lon_lat[:,1] = particles_lat_all[time_dex,:]
        else:
            points_lon_lat[:,0] = particles_lon_all[:,time_dex]
            points_lon_lat[:,1] = particles_lat_all[:,time_dex]

        start_1 = time.time()

        for polygon_dex in range(num_polygons):

            path = plt_path.Path(polygon_vertex_list[polygon_dex])
            
            particles_inside_flags = path.contains_points(points_lon_lat)
       
            settler_indices = np.logical_and(particles_inside_flags,particle_safety_mask)  

            settlement_boxes[settler_indices] = polygon_dex + 1
            settlement_times[settler_indices] = time_dex + 1
            
            # Update the counter of settled particles.  For checking consistency
            counter_array[:] += settler_indices 
            #counter_array[:,file_number] += settler_indices 

            # Update the safety mask, so that settlement prevents further modifications to the stored settlement location
            particle_safety_mask[particles_inside_flags] = False
            
        
        end_1 = time.time()
        print(f'slow loop took {(end_1-start_1)/60} minutes')

    start = time.time()

    
    valid_timesteps = particles_lon_all != dummy_value # For time averages

    if time_dimension == 0:
        particle_drift_times = np.sum(valid_timesteps[0:timesteps_full_run,:], axis=1)
    else:
        particle_drift_times = np.sum(valid_timesteps[:,0:timesteps_full_run], axis=1)
        
    if time_dimension == 1:
        #particle_arrays_O2 = np.sum(O2_threshold_breach_array[:,:,0:timesteps_full_run], axis=2)
        #particle_arrays_pH = np.sum(pH_threshold_breach_array[:,:,0:timesteps_full_run], axis=2)
       
        # Compute average temperature along trajectories
        particle_array_T = np.sum(particles_T_all[:,0:timesteps_full_run], axis=1) 
        average_T = np.divide(particle_array_T,particle_drift_times)
        average_T[np.isinf(average_T)] = np.nan   # and just bin this later with hist.  make sure to use plt.hist(average_T[np.logical_not(np.isnan(average_T))]
        average_T_round = np.round(average_T).astype(np.int64)
        average_T_round = np.where(np.abs(average_T_round) > dummy_value, np.nan, average_T_round)
        average_T_scaled = average_T_round * T_scale_factor
   
    #print('ending slow loop')
    
    end = time.time()

    print(f'fast loop took {(end-start)/60} minutes')




    num_good = 0
    num_bad = 0
    num_empty = 0
    num_outofbounds = 0

    bad_release_count = 0

    # modify the pdf data structure
    for particle_num in range(num_particles):

        if release_boxes[particle_num] < 0:
            raise RuntimeError('particle {} has release cell dex: {}'.format(particle_num,release_boxes[particle_num]))
        # --------------------------------------------------
        # --------------------------------------------------
            
        # Patrick only has the annual array
        array_dex_list = [0]
       
        
        # wait we have 0's in here, so my logic requires we skip in this case
        if release_boxes[particle_num] == 0:
            if time_dimension == 0:
                if particles_lon_all[0,particle_num] == dummy_value:
                    num_empty += 1
            else:
                if particles_lon_all[particle_num,0] == dummy_value:
                    num_empty += 1
            num_bad += 1
            continue
        #elif good_indices[particle_num]:
        else:
            num_good += 1
            release_counts_per_cell[0,int(release_boxes[particle_num])-1] += 1
            for array_dex in array_dex_list:

                if settlement_boxes[particle_num] > 0:

                    pdf_arrays_connectivity[array_dex,int(release_boxes[particle_num])-1,int(settlement_boxes[particle_num])-1] += 1   
           
                    if time_dimension == 1:

                        pdf_arrays_T[array_dex,int(release_boxes[particle_num])-1,int(average_T_scaled[particle_num])] += 1

                        #for jj in range(len(O2_limit_list)):
                        #    if particle_arrays_O2[jj,particle_num] > 0:
                        #        pdf_arrays_O2[jj,array_dex,int(release_boxes[particle_num])-1,int(particle_arrays_O2[jj,particle_num])-1] += 1

                        #for jj in range(len(pH_limit_list)):
                        #    if particle_arrays_pH[jj,particle_num] > 0:
                        #        pdf_arrays_pH[jj,array_dex,int(release_boxes[particle_num])-1,int(particle_arrays_pH[jj,particle_num])-1] += 1

        

    print('num_good: {}'.format(num_good))
    print('num_empty: {}'.format(num_empty))
    print('num_outofbounds: {}'.format(num_outofbounds))
    print('num_bad: {}'.format(num_bad))

    # Save the stationary settlers counts!
    stationary_settlers_array = stationary_settlers_per_box
    #stationary_settlers_array[:] = stationary_settlers_per_box
    #stationary_settlers_array[:,file_number] = stationary_settlers_per_box
    # Save the number of empty trajectories in a file
    #empty_trajectory_count_array = num_empty
    empty_trajectory_count_array[file_number] = num_empty

    file_number += 1

    end_file_time = time.time()
    print(f'File processing took {(end_file_time-start_file_time)/60} minutes')
   
    # End of loop that used to be here
    # ------------------------------------------------------------------------------------------------
    #

    # Remove stationary local settlers
    num_stationary_settlers = 0
    for polygon_num in range(num_polygons):
        num_stationary_settlers_polygon = int(np.sum(stationary_settlers_array[polygon_num,:]))
        #print(f'Number of stationary settlers: {num_stationary_settlers_polygon}')
        pdf_arrays_connectivity[array_dex,polygon_num,polygon_num] -= num_stationary_settlers_polygon
        num_stationary_settlers += num_stationary_settlers_polygon
        #pdf_arrays_connectivity[array_dex,polygon_num,polygon_num] -= int(stationary_settlers_per_box[polygon_num])

    print('Number of stationary settlers: {}'.format(num_stationary_settlers))

    d = {}
    d['release_counts_per_cell'] = release_counts_per_cell
    d['pdf_arrays_connectivity'] = pdf_arrays_connectivity
    d['pdf_arrays_settleTime'] = pdf_arrays_settleTime
    d['counter_array'] = counter_array
    d['stationary_settlers_array'] = stationary_settlers_array
    d['empty_trajectory_count_array'] = empty_trajectory_count_array
    d['pld_days'] = pld_array[pld_dex]
    
    if time_dimension == 1:
        d['pdf_arrays_T'] = pdf_arrays_T
        #d['pdf_arrays_O2'] = pdf_arrays_O2
        #d['pdf_arrays_pH'] = pdf_arrays_pH
        #d['O2_limit_list'] = O2_limit_list
        #d['pH_limit_list'] = pH_limit_list

    save_output_file_name_pre = "conn_hist_data_"
    save_output_file_name = save_output_file_name_pre + tracking_file.split('.')[0].split('/')[-1]
    #save_output_file_name = save_output_file_name_pre + tracking_file.split('.')[0]
    save_output_full_path = os.path.join(save_output_directory, save_output_file_name + f"_pld_{pld_array[pld_dex,0]}_{pld_array[pld_dex,1]}_version_{script_version}_{polygon_string}")
    

    np.savez(save_output_full_path, **d)







