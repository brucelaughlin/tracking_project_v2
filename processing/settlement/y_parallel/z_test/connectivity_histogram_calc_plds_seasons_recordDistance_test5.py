# NOTE TO USER:
# "prePdf" in a variable name indicates that the variable will later be normalized; I didn't want to call these variables "pdfs" because
# they aren't pdfs yet, and I didn't know of a succint name for the "gridded data" structures that I'm creating which will become pdfs.  Hence, prePdfs!

#script_version = "basicVars"
script_version = "recordDistances_test5"

# testing with this dir:
# /data03/blaughli/tracking_output/z_complete/z_tests/Mercator_coastalCells_1993_1993_kick_0p05_singleRun_10days___Mercator



# Opendrift output files have the following for time units:
# time:units = "seconds since 1970-01-01"
#base_year = 1970


# Hardcoding these paths, JUST FOR TESTING  
#tracking_file = "/home/blaughli/tracking_project_v2/t_scraps/test_dir_con/tracking_output_configFile_000_job_00.nc"
##tracking_file = "/home/blaughli/tracking_project_v2/t_scraps/dummy_test_dir_con/a_test_file.nc"
##tracking_file = "/home/blaughli/tracking_project_v2/t_scraps/dummy_test_dir_con_10/a_test_file_10.nc"

#polygon_file_path = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/bounding_boxes_lonlat_Mercator_singleCoastalCells.txt"



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
from geopy.distance import geodesic

start_script_time = time.time()

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# PARAMETERS and CONSTANTS

parser = argparse.ArgumentParser()
parser.add_argument("--trackingfile", type=str)
parser.add_argument("--pdrakefileswitch", type=int)
parser.add_argument("--polygonfile", type=str)
#parser.add_argument("--baseyear", type=str)
args = parser.parse_args()

#tracking_file = "/data03/blaughli/tracking_output/Mercator_coastalCells_1993_2018_kickSTD_0p0___global-reanalysis-phy-001-030-daily_1993_2018/tracking_output_configFile_000_job_00.nc"
#pdrake_file_switch = 0
#polygon_file_path = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/bounding_boxes_lonlat_Mercator_singleCoastalCells.txt"



tracking_file = args.trackingfile
pdrake_file_switch = args.pdrakefileswitch
polygon_file_path = args.polygonfile


tracking_dir = os.path.dirname(tracking_file)

polygon_string = Path(polygon_file_path).stem


# Opendrift swaps the particle and time dimensions relative to the tracking code that Patrick used
# For my files, the particle dimension is the first dimension (ie 0):
particle_dimension = 0
time_dimension = 1
# For Patrick's files, the particle dimension is the second dimension (ie 1):
if pdrake_file_switch != 0:
    particle_dimension = 1
    time_dimension = 0



base_datetime = datetime.datetime(base_year,1,1,0,0,0)



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

#pld_list = np.array([[5,6],[10,11],[15,17],[20,22],[30,33],[45,49],[60,65],[90,98],[120,131]])
#pld_list=np.array([pld_kelp_bass,pld_ca_sheephead,pld_kelp_rockfish,pld_blue_black_rockfish,pld_test])
#pld_list=np.array([pld_kelp_bass,pld_ca_sheephead,pld_kelp_rockfish,pld_blue_black_rockfish])

#pld_list = [[5,6],[10,11],[15,17],[20,22],[30,33],[45,49],[60,65],[90,98],[120,131]]

# Pete Raimondi's PLDs
pld_list = [[1,5],[5,10],[10,15],[15,30],[30,45],[45,60],[60,75],[75,90],[90,120],[120,150],[150,180]]

# This is for my clunky idea of first determining all polygons all floats enter
last_pld_day_global = pld_list[-1][1]


# choose one pld to use, while still testing
#pld_chosen_dex = num_plds-1
#pld_chosen_dex = 2
pld_chosen_dex = 5


num_plds = len(pld_list)

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

save_output_directory = tracking_dir + 'v_connHist_files_test5'
#save_output_directory = tracking_dir + 'v_connHist_files'

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

# number of prePdfs per feild (4 seasonal and 1 annual = 5)
# JUST MAKE THE ANNUAL FILE
#num_prePdfs = 1
# JUST KIDDING 
num_prePdfs = 5

# Store the index of the annual prePdf used in the connectivity array (original idea is seasons are indices 0-3, and annual is index 4)
annual_prePdf_index = num_prePdfs - 1


# ----------------------------------------------------------loading
dset = netCDF4.Dataset(tracking_file, 'r')

if time_dimension == 0:
    ocean_time = dset.variables['ocean_time'][:]
    status_all = np.zeros_like(lon_all.T)
else:
    ocean_time = dset.variables['time'][:]
    status_all = dset.variables['status'][:]
    # Exposure variables 
    #O2_all = dset.variables['oxygen'][:]
    #pH_all = dset.variables['pH'][:]
#    temp_all = dset.variables['sea_water_temperature'][:]
    
    #O2_all *= conversion_factor
#counter_array=np.zeros((num_particles))

lon_all = dset.variables['lon'][:]
lat_all = dset.variables['lat'][:]


dset.close()
# ----------------------------------------------------------loading


num_particles = np.shape(lon_all)[particle_dimension]


print('NUMBER OF PARTICLES: {}'.format(num_particles))


timesteps_per_day = int(1/((ocean_time[1]-ocean_time[0])/86400))
#timesteps_per_day = int(1/((ocean_ocean_time[1]-ocean_ocean_time[0])/86400))

timesteps_full_run = len(ocean_time) * timesteps_per_day

# Here again assuming we'll never handle timesteps bigger than 1 day (see the "1" below)
timestep_length = 1/timesteps_per_day

# For now, exit if output is less than daily...
if timesteps_per_day < 1:
    raise RuntimeError('This code is dumb and currently breaks if the save timestep is > 1 day.  Should be an easy fix')



## NEW ALGORITHM
pld_list_timesteps = pld_list.copy()
for pld_dex in range(len(pld_list)):
    for day_dex in range(2):
        pld_list_timesteps[pld_dex][day_dex] *= timesteps_per_day
        #print(pld_list_timesteps[pld_dex][day_dex])



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
# Create lists to store statistics for full run and seasonal subsets
#---------------------------------------------------------------------
#### THIS USED TO BE INSIDE THE LOOP... BUT NOW IT'S OUTSIDE, SINCE WE'RE PROCESSING JUST A SINGLE FILE IN THIS SCRIPT

# store the number of particles released from each box
release_counts_per_polygon_array = np.zeros((num_prePdfs,num_polygons), dtype=int)
#release_counts_per_polygon_array = np.zeros((num_prePdfs,num_plds,num_polygons), dtype=int)

# Connectivity (release box number vs settlement box number)
prePdf_arrays_connectivity = np.zeros((num_prePdfs,num_plds,num_polygons,num_polygons))
prePdf_arrays_connectivity_noStagnation = np.zeros((num_prePdfs,num_plds,num_polygons,num_polygons))

particle_distances_per_pld_settlers = np.zeros((num_particles,num_plds))
particle_distances_per_pld_settlers_noStagnation = np.zeros((num_particles,num_plds))
#particle_distances_per_pld_settlers = np.full((num_particles,num_plds), np.nan)
#particle_distances_per_pld_settlers_noStagnation = np.full((num_particles,num_plds), np.nan)
particle_distances_per_pld_allTimesteps = np.zeros((num_particles,num_plds))

polygons_settled_per_particle = np.zeros((num_particles,num_plds))
stagnant_settlers_per_particle_per_pld = np.zeros((num_particles,num_plds))

particle_settle_times_per_pld = np.zeros((num_particles,num_plds))
particle_settle_times_per_pld_noStagnation = np.zeros((num_particles,num_plds))

## Time after PLD until settlement (saving only release location) (release box number vs settlement time)
#prePdf_arrays_settleTime = np.zeros((num_prePdfs,num_plds,num_polygons,timesteps_settlement_window))
#prePdf_arrays_settleTime_noStagnation = np.zeros((num_prePdfs,num_plds,num_polygons,timesteps_settlement_window))

## Number of days eposed to DO levels below 2.2 (saving only settlement location) (release box number vs exposure time)
## note: Think the time dimension needs to be one bigger than the number of possible timesteps (ie need to include an option for "zero") 
#prePdf_arrays_O2 = np.zeros((len(O2_limit_list),num_prePdfs,num_plds,num_polygons,timesteps_full_run))

## Same idea for pH
#prePdf_arrays_pH = np.zeros((len(pH_limit_list),num_prePdfs,num_plds,num_polygons,timesteps_full_run))

## Histogram of average temperature experienced - just estimating range, since I don't know it without processing
#prePdf_arrays_T = np.zeros((num_prePdfs,num_plds,num_polygons,n_T_bins))



stationary_settlers_array = np.zeros((num_polygons,1))

empty_trajectory_count_array = np.zeros(1)

#print(tracking_file)



# Store the total number of particles
num_particles = np.shape(lon_all)[particle_dimension]
num_timesteps = np.shape(lon_all)[time_dimension]

# SEASONS - Prepare the list of possible seed months for the run
run_seed_months_list = []
for t in ocean_time:
    run_seed_months_list.append(datetime.datetime.strptime(str(base_datetime+datetime.timedelta(seconds=t)), '%Y-%m-%d %H:%M:%S').month)

# Preapre list for storing the months in which seeding happened...?
#seed_months = []
# Update: use an array, and initialize to -1 so that -1 can be ignored later
seed_months_per_particle_array = np.ones(num_particles) * -1

# Why not just make one array of all of the data, and then reference it by time index???
if time_dimension == 0:
    particles_lon_all=np.full((timesteps_full_run,num_particles), np.nan)
    particles_lat_all=np.full((timesteps_full_run,num_particles), np.nan)
    #particles_O2_all=np.zeros([timesteps_full_run,num_particles])
    #particles_pH_all=np.zeros([timesteps_full_run,num_particles])
#    particles_T_all=np.zeros([timesteps_full_run,num_particles])
else: 
    particles_lon_all=np.full((num_particles,timesteps_full_run), np.nan)
    particles_lat_all=np.full((num_particles,timesteps_full_run), np.nan)
    #particles_O2_all=np.zeros([num_particles,timesteps_full_run])
    #particles_pH_all=np.zeros([num_particles,timesteps_full_run])
#    particles_T_all=np.zeros([num_particles,timesteps_full_run])

num_invalid_status = 0


for particle_num in range(num_particles):

    # Need trajectory status (active vs inactive) for determining month of release of particles
    trajectory_status = status_all[particle_num,:]
    #trajectory_mask = trajectory_status == 0

    # Store month of first timestep of trajectory (ie seeding/starting month)

    # NOTE: Had cases where a particle had inactive status (1) for its whole life... not sure why, but I had clearly assumed this wouldn't
    # happen (it seems sketchy now to be using a list, appending a value for each particle.  why not just start with a fixed array, and
    # for problematic particles, leave their value as the default, ie nan?

    if np.shape(np.where(trajectory_status == 0))[1] == 0:
        num_invalid_status += 1
        continue
    #seed_months.append(run_seed_months_list[np.where(trajectory_status == 0)[0][0]])
    seed_months_per_particle_array[particle_num] = run_seed_months_list[np.where(trajectory_status == 0)[0][0]]

    if time_dimension == 0:
        particle_lon = lon_all[:,particle_num]
        particle_lon = particle_lon[np.logical_not(particle_lon.mask)].data
        #particle_lon = particle_lon[0:timesteps_full_run]
        # ADDING A STEP WHICH IS NECESSARY WHEN THERE IS NO MASKED DATA (CALLING FOR THE DATA ABOVE ADDS A DIMENSION, WHICH MESSES UP INDEXING)
        if np.shape(np.shape(particle_lon))[0] == 2:
            particle_lon = particle_lon[0,:]
        particle_lat = lat_all[:,particle_num]
        particle_lat = particle_lat[np.logical_not(particle_lat.mask)].data
        #particle_lat = particle_lat[0:timesteps_full_run]
        if np.shape(np.shape(particle_lat))[0] == 2:
            particle_lat = particle_lat[0,:]
    else:
        particle_lon = lon_all[particle_num,:]
        particle_lon = particle_lon[np.logical_not(particle_lon.mask)].data
        #particle_lon = particle_lon[0:timesteps_full_run]
        if np.shape(np.shape(particle_lon))[0] == 2:
            particle_lon = particle_lon[0,:]
        particle_lat = lat_all[particle_num,:]
        particle_lat = particle_lat[np.logical_not(particle_lat.mask)].data
        #particle_lat = particle_lat[0:timesteps_full_run]
        if np.shape(np.shape(particle_lat))[0] == 2:
            particle_lat = particle_lat[0,:]
    
        #particle_O2 = O2_all[particle_num,:]
        #particle_O2 = particle_O2[np.logical_not(particle_O2.mask)].data
        #particle_O2 = particle_O2[0:timesteps_full_run]
        #particle_pH = pH_all[particle_num,:]
        #particle_pH = particle_pH[np.logical_not(particle_pH.mask)].data
        #particle_pH = particle_pH[0:timesteps_full_run]
 #       particle_T = temp_all[particle_num,:]
 #       particle_T = particle_T[np.logical_not(particle_T.mask)].data
 #       particle_T = particle_T[0:timesteps_full_run]
    
    if time_dimension == 0:
        particles_lon_all[0:len(particle_lon),particle_num] = particle_lon
        particles_lat_all[0:len(particle_lat),particle_num] = particle_lat
    else:
        particles_lon_all[particle_num,0:len(particle_lon)] = particle_lon
        particles_lat_all[particle_num,0:len(particle_lat)] = particle_lat
        #particles_O2_all[particle_num,:] = np.pad(particle_O2,(0,timesteps_full_run-len(particle_O2)), 'constant',constant_values=(dummy_value))
        #particles_pH_all[particle_num,:] = np.pad(particle_pH,(0,timesteps_full_run-len(particle_pH)), 'constant',constant_values=(dummy_value))
#        particles_T_all[particle_num,:] = np.pad(particle_T,(0,timesteps_full_run-len(particle_T)), 'constant',constant_values=(dummy_value))

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
    
#    particles_T_all = np.where(particles_T_all == dummy_value, np.nan, particles_T_all)
#    particles_T_all = np.ma.array(particles_T_all, mask = np.isnan(particles_T_all))
    
    #O2_threshold_breach_array = np.where(O2_threshold_breach_array == dummy_value, np.nan, O2_threshold_breach_array)
    #O2_threshold_breach_array = np.ma.array(O2_threshold_breach_array, mask = np.isnan(O2_threshold_breach_array))
    
    #pH_threshold_breach_array = np.where(pH_threshold_breach_array == dummy_value, np.nan, pH_threshold_breach_array)
    #pH_threshold_breach_array = np.ma.array(pH_threshold_breach_array, mask = np.isnan(pH_threshold_breach_array))
    
    # Exposure statistics
    #particle_arrays_O2 = np.zeros((num_particles,len(O2_limit_list)), dtype=int)
    #particle_arrays_pH = np.zeros((num_particles,len(pH_limit_list)), dtype=int)
#    particle_array_T = np.zeros(num_particles)
#    particle_array_driftTime = np.zeros(num_particles)



# --------------------------------------------------------------------
# WHY WAS THIS REDUNDANT CODE HERE?
# --------------------------------------------------------------------
# LOOK AT NP.MAWHERE - ACCOMPLISH IN ONE STEP
# or np.maand - may be an "and" for masked arrays, which takes masks into account

#particles_lon_all = np.where(particles_lon_all == dummy_value, np.nan, particles_lon_all)
#particles_lon_all = np.ma.array(particles_lon_all, mask = np.isnan(particles_lon_all))

#particles_lat_all = np.where(particles_lat_all == dummy_value, np.nan, particles_lat_all)
#particles_lat_all = np.ma.array(particles_lat_all, mask = np.isnan(particles_lat_all))
# --------------------------------------------------------------------


# We need to track a particle's current polygon location at all times, to remove particles which never leave their seeding polygon
# New idea: just have a single boolean vector of length <num_particles> which starts as true (ie assumes particle has not left seed polygon),
# and switches to false when the particle leaves the seed polygon, and never changes value again.  So, we start by assuming they're all stagnant
stagnant_particle_boolean_vector = np.ones(num_particles, dtype=bool)
#stagnant_particle_boolean_vector = np.zeros(num_particles, dtype=bool)
#stagnant_particle_boolean_vector = np.zeros(num_particles)==0
#particles_polygons_full_run = np.zeros_like(particles_lon_all)
#particles_polygons_full_run[:] = np.nan

# Done later
#### Create array for tracking drift time, for use in computing average T along trajectory 
###particle_drift_times = np.zeros(num_particles)

# prepare lists to hold the starting/ending polygon numbers for each particle
release_polygons = np.zeros(num_particles, dtype=int)

#settlement_polygons = np.zeros(num_particles, dtype=int)
#settlement_polygons_noStagnation = np.zeros(num_particles, dtype=int)
## add settlement time storage
#settlement_times = np.zeros(num_particles, dtype=int)
#settlement_times_noStagnation = np.zeros(num_particles, dtype=int)

settlement_polygons = np.zeros((num_particles,num_plds), dtype=int)
settlement_polygons_noStagnation = np.zeros((num_particles,num_plds), dtype=int)
settlement_times = np.zeros((num_particles,num_plds), dtype=int)
settlement_times_noStagnation = np.zeros((num_particles,num_plds), dtype=int)


#---------------------------------------------------------------------
# FIRST MUST DETERMINE STARTING LOCATIONS
#---------------------------------------------------------------------
unsettled_particle_mask = np.ones(num_particles, dtype=bool)

# <points_lon_lat> is used here and in settlement; for starting locations, its value is fixed

points_lon_lat = np.zeros((num_particles,2))
if time_dimension == 0:
    points_lon_lat[:,0] = particles_lon_all[0,:]
    points_lon_lat[:,1] = particles_lat_all[0,:]
else:
    points_lon_lat[:,0] = particles_lon_all[:,0]
    points_lon_lat[:,1] = particles_lat_all[:,0]

# Actually, beyond the counts of releases per polygon for use in making pdfs later, I now also want a vector with the starting polygon (with my +1 convention,
# so the first is number 1 (not 0)) for each particle.  I'll try to use this for determining when particles stop being stagnant
particle_seed_polygon_vector = np.zeros(np.shape(points_lon_lat)[0])


for polygon_dex in range(num_polygons):
    path = plt_path.Path(polygon_vertex_list[polygon_dex])
    particles_inside_flags = path.contains_points(points_lon_lat)

    particle_seed_polygon_vector[particles_inside_flags] = polygon_dex + 1

    release_polygons[np.logical_and(particles_inside_flags,unsettled_particle_mask)] = polygon_dex + 1
    unsettled_particle_mask[particles_inside_flags] = False


# Need a mask to prevent processing of particles with no release location (I have no idea how that happens, by the way..???!)
empty_release_cells = np.setdiff1d(np.arange(num_polygons)+1,release_polygons)

print("--------------------------------------------")
print(f"{len(empty_release_cells)} Empty release cells:")
for ii in range(len(empty_release_cells)):
    print(f"cell {empty_release_cells[ii]}")
print("--------------------------------------------")


#particles_polygons_full_run[0,:] = release_polygons

# -------------------------------------------------------------------------------------------------------------------
# Starting loop below this.  Old version was already in loop, so the code above was already indented
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------

    
# NOW DETERMINE SETTLEMENT LOCATIONS!
# NOTE: Release locations could probably be easily calculated as part of the main loop (ie timestep 0), but for now I don't
# want to mess anything up so I'm leaving that as a separate calculation.  



#---------------------------------------------------------------------
# START THE MAIN LOOP!!!
#---------------------------------------------------------------------

points_lon_lat = np.zeros((num_particles,2))
points_lon_lat_previous = np.zeros((num_particles,2))

# We're now gonna determine the pld dex dynamically, as we loop over ALL output timesteps
pld_dex = 0

currently_within_pld_switch = False

first_time_dex = 0

first_timestep_of_run_switch = True

# Create the "safety mask" that I'll use to make sure I only store the first polygon entered during the settlement window
# (not used if current timestep is not within a pld)
unsettled_particle_mask = np.ones(num_particles, dtype=bool)

for time_dex in range(first_time_dex,pld_list_timesteps[-1][-1]):

    print('timestep {}/{}'.format(time_dex,pld_list_timesteps[-1][-1] - 1))

    #first_timestep_of_run_switch = time_dex == first_time_dex

    last_timestep_of_pld_switch = False
    
    # create the "safety mask" that I'll use to make sure I only store the first polygon entered during the settlement window
    # (not used if current timestep is not within a pld)
#    unsettled_particle_mask = np.ones(num_particles, dtype=bool)
    
    if time_dex >= pld_list_timesteps[pld_dex][0] and time_dex < pld_list_timesteps[pld_dex][1]:
        currently_within_pld_switch = True
        # Need switch for indicating last timestep of pld, for storing connectivity
        if time_dex == pld_list_timesteps[pld_dex][1] - 1:
            last_timestep_of_pld_switch = True

    elif time_dex >= pld_list_timesteps[pld_dex][1]:
        unsettled_particle_mask = np.ones(num_particles, dtype=bool) # I hadn't been re-setting the safety mask with each new PLD... ???!?!?!?!?!?!?
        pld_dex += 1
        if pld_dex == num_plds:
            break
        elif time_dex >= pld_list_timesteps[pld_dex][0]:
            currently_within_pld_switch = True
        else:
            currently_within_pld_switch = False
    #else:
    #    currently_within_pld_switch = False
    
#    print(f"Currently within PLD: {currently_within_pld_switch}")

    '''
    # to handle a pld, eventually want something like:
    if time_dex >= pld_list_timesteps[pld_dex][0] and time_dex < pld_list_timesteps[pld_dex][1]:
        new_pld_switch = True
        ...settlement stuff

    elif new_pld_switch:
        pld_dex += 1
        new_pld_switch = False
    '''
    
    #for pld_dex in range(num_plds):
        
    ###TESTING
    # test using the chosen pld
    #if pld_dex != pld_chosen_dex:
    #    continue
        #break


#    # Prepare array for counting settlers which never leave their seed polygon
#    stationary_settlers_per_polygon = np.zeros((num_polygons,1))
   
#    first_settlement_day = pld_list[pld_dex,0]
#    last_settlement_day = pld_list[pld_dex,1]

#    print('pld: {}-{}'.format(first_settlement_day,last_settlement_day))

    # THIS is the v5 adjustment: need to only use data from the pld, which begins after "first_settlement_day" and ends after "last_settlement_day"
#    pld_length_days = last_settlement_day - first_settlement_day + 1


#    timesteps_settlement_window = pld_length_days * timesteps_per_day
#    timesteps_full_run = (last_settlement_day+1) * timesteps_per_day + 1 # why isn't this just taken from the dimensions of the files?


    #first_settle_dex = first_settlement_day * timesteps_per_day + 1
#    first_settle_dex = first_settlement_day * timesteps_per_day  # Why was I adding 1 here????
#    last_settle_dex = (last_settlement_day+1) * timesteps_per_day # And here????
    #last_settle_dex = (last_settlement_day+1) * timesteps_per_day + 1



#    stationary_settler_array = np.ones((num_particles,num_polygons), dtype=bool) 

    #for time_dex in range(first_settle_dex,last_settle_dex):

    #active_current_polygon_occupation = np.empty(num_particles)
    #active_current_polygon_occupation[:] = np.nan
    active_current_polygon_occupation = np.zeros(num_particles)


    if time_dex > first_time_dex:
        points_lon_lat_previous = points_lon_lat.copy()

    if time_dimension == 0:
        points_lon_lat[:,0] = particles_lon_all[time_dex,:]
        points_lon_lat[:,1] = particles_lat_all[time_dex,:]
    else:
        points_lon_lat[:,0] = particles_lon_all[:,time_dex]
        points_lon_lat[:,1] = particles_lat_all[:,time_dex]

    start_1 = time.time()

    # must reverse coordinates, since lat/lon is standard but I've been using lon/lat
    if time_dex > first_time_dex:
        for particle_dex in range(num_particles):
            
            point1 = list(points_lon_lat_previous[particle_dex])
            point1.reverse()
            point2 = list(points_lon_lat[particle_dex])
            point2.reverse()
           
            if np.sum(np.isnan(point1)) > 0 or np.sum(np.isnan(point2)):
                continue

            distance = geodesic(point1, point2).km  # Returns distance in kilometers
            particle_distances_per_pld_allTimesteps[particle_dex,pld_dex:] += distance
            #particle_distances_per_pld_allTimesteps[particle_dex,pld_dex] += distance

    
#    if time_dex > 12:
#        break



    for polygon_dex in range(num_polygons):

        path = plt_path.Path(polygon_vertex_list[polygon_dex])
        
        particles_inside_flags = path.contains_points(points_lon_lat)
  
        stagnant_particle_boolean_vector[particles_inside_flags] = np.logical_and(particle_seed_polygon_vector[particles_inside_flags] == polygon_dex + 1, stagnant_particle_boolean_vector[particles_inside_flags])
        #stagnant_particle_boolean_vector[particles_inside_flags] = np.logical_and(particle_seed_polygon_vector[particles_inside_flags] == polygon_dex + 1, np.logical_not(stagnant_particle_boolean_vector[particles_inside_flags]))

        active_current_polygon_occupation_indices = np.logical_and(particles_inside_flags, unsettled_particle_mask)  
        active_current_polygon_occupation_indices_noStagnation = np.logical_and(np.logical_and(particles_inside_flags, unsettled_particle_mask), np.logical_not(stagnant_particle_boolean_vector))  
        #active_current_polygon_occupation_indices_noStagnation = np.logical_and(np.logical_and(particles_inside_flags, unsettled_particle_mask), stagnant_particle_boolean_vector)  

        active_current_polygon_occupation[active_current_polygon_occupation_indices] = polygon_dex + 1

        if currently_within_pld_switch:

            particle_distances_per_pld_settlers[active_current_polygon_occupation_indices,pld_dex] = particle_distances_per_pld_allTimesteps[active_current_polygon_occupation_indices,pld_dex]
            particle_distances_per_pld_settlers_noStagnation[active_current_polygon_occupation_indices_noStagnation,pld_dex] = (
                    particle_distances_per_pld_allTimesteps[active_current_polygon_occupation_indices_noStagnation,pld_dex])
            
            polygons_settled_per_particle[active_current_polygon_occupation_indices,pld_dex] = polygon_dex + 1

            #stagnant_settlers_per_particle_per_pld[:,pld_dex] += stagnant_particle_boolean_vector
            stagnant_settlers_per_particle_per_pld[active_current_polygon_occupation_indices,pld_dex] += stagnant_particle_boolean_vector[active_current_polygon_occupation_indices].astype(int)

            particle_settle_times_per_pld[active_current_polygon_occupation_indices,pld_dex] = time_dex * timestep_length 
#            particle_settle_times_per_pld_noStagnation



            settlement_polygons[active_current_polygon_occupation_indices,pld_dex] = polygon_dex + 1
            settlement_times[active_current_polygon_occupation_indices,pld_dex] = time_dex + 1

            settlement_polygons_noStagnation[active_current_polygon_occupation_indices_noStagnation,pld_dex] = polygon_dex + 1
            settlement_times_noStagnation[active_current_polygon_occupation_indices_noStagnation,pld_dex] = time_dex + 1
            
            # Update the safety mask, so that settlement prevents further modifications to the stored settlement location
            unsettled_particle_mask[particles_inside_flags] = False

#        particles_polygons_full_run[time_dex,:] = settlement_polygons
        


        # Update the counter of settled particles.  For checking consistency
#        counter_array[:] += active_current_polygon_occupation_indices 
        #counter_array[:,file_number] += active_current_polygon_occupation_indices 


#        stagnant_particle_boolean_vector = np.logical_and(release_polygons == active_current_polygon_occupation, stagnant_particle_boolean_vector)
        #stagnant_particle_boolean_vector = np.logical_and(np.logical_and(release_polygons, active_current_polygon_occupation), stagnant_particle_boolean_vector)

    #print(f"Number of stagnant particles: {np.sum(stagnant_particle_boolean_vector)}")

    
    end_1 = time.time()
#    print(f'slow loop took {(end_1-start_1)/60} minutes')
    

    











    
#    start = time.time()


#    valid_timesteps = particles_lon_all != dummy_value # For time averages
#
#    if time_dimension == 0:
#        particle_drift_times = np.sum(valid_timesteps[0:timesteps_full_run,:], axis=0) #right? axis 0 not axis 1?
#        #particle_drift_times = np.sum(valid_timesteps[0:timesteps_full_run,:], axis=1)
#    else:
#        particle_drift_times = np.sum(valid_timesteps[:,0:timesteps_full_run], axis=1)
#        
#    if time_dimension == 1:
#        #particle_arrays_O2 = np.sum(O2_threshold_breach_array[:,:,0:timesteps_full_run], axis=2)
#        #particle_arrays_pH = np.sum(pH_threshold_breach_array[:,:,0:timesteps_full_run], axis=2)
#       
#        # Compute average temperature along trajectories
#        particle_array_T = np.sum(particles_T_all[:,0:timesteps_full_run], axis=1) 
#        average_T = np.divide(particle_array_T,particle_drift_times)
#        average_T[np.isinf(average_T)] = np.nan   # and just bin this later with hist.  make sure to use plt.hist(average_T[np.logical_not(np.isnan(average_T))]
#        average_T_round = np.round(average_T).astype(np.int64)
#        average_T_round = np.where(np.abs(average_T_round) > dummy_value, np.nan, average_T_round)
#        average_T_scaled = average_T_round * T_scale_factor
   
    #print('ending slow loop')
    
#    end = time.time()

#    print(f'fast loop took {(end-start)/60} minutes')


#    # Create the current pld's connectivity prePdf if current timestep is the last timestep of the current pld
#    if last_timestep_of_pld_switch:
                    
                    
                    
#firstPldSwitch = False



firstPldSwitch = True
for pld_dex in range(num_plds):

    num_good = 0
    num_bad = 0
    num_empty = 0

    bad_release_count = 0

    # modify the prePdf data structure
    for particle_num in range(num_particles):

        # --------------------------------------------------
        if release_polygons[particle_num] < 0:
            raise RuntimeError('particle {} has release cell dex: {}'.format(particle_num,release_polygons[particle_num]))
        # --------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Seasons
        # We always add to the "annual" array; and we also add to the seasonal array
        # -------------------------------------------------------------------------------------
        if seed_months_per_particle_array[particle_num] >=0 and seed_months_per_particle_array[particle_num] <=2 :
            season_dex = 0 #DJF
            #season_dex = 1 #DJF
        elif seed_months_per_particle_array[particle_num] >=3 and seed_months_per_particle_array[particle_num] <=5 :
            season_dex = 1 #MAM
            #season_dex = 2 #MAM
        elif seed_months_per_particle_array[particle_num] >=6 and seed_months_per_particle_array[particle_num] <=8 :
            season_dex = 2 #JJA
            #season_dex = 3 #JJA
        elif seed_months_per_particle_array[particle_num] >=9 and seed_months_per_particle_array[particle_num] <=11 :
            season_dex = 3 #SON
            #season_dex = 4 #SON
        #else:

        array_dex_list = [season_dex,annual_prePdf_index]
        # -------------------------------------------------------------------------------------


        #### Patrick only has the annual array
        #array_dex_list = [0]
       
        
        # wait we have 0's in here, so my logic requires we skip in this case
        if release_polygons[particle_num] == 0:
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
            #release_counts_per_polygon_array[0,int(release_polygons[particle_num])-1] += 1
            for array_dex in array_dex_list:
                if firstPldSwitch:
                    release_counts_per_polygon_array[array_dex,int(release_polygons[particle_num])-1] += 1
                if settlement_polygons[particle_num,pld_dex] > 0:
                    prePdf_arrays_connectivity[array_dex,pld_dex,int(release_polygons[particle_num])-1,int(settlement_polygons[particle_num,pld_dex])-1] += 1   
                if settlement_polygons_noStagnation[particle_num,pld_dex] > 0:
                    prePdf_arrays_connectivity_noStagnation[array_dex,pld_dex,int(release_polygons[particle_num])-1,int(settlement_polygons[particle_num,pld_dex])-1] += 1   
                
                
                    #if time_dimension == 1:

                    #    prePdf_arrays_T[array_dex,pld_dex,int(release_polygons[particle_num])-1,int(average_T_scaled[particle_num])] += 1
                    #    #for jj in range(len(O2_limit_list)):
                    #    #    if particle_arrays_O2[jj,particle_num] > 0:
                    #    #        prePdf_arrays_O2[jj,array_dex,int(release_polygons[particle_num])-1,int(particle_arrays_O2[jj,particle_num])-1] += 1

                    #    #for jj in range(len(pH_limit_list)):
                    #    #    if particle_arrays_pH[jj,particle_num] > 0:
                    #    #        prePdf_arrays_pH[jj,array_dex,int(release_polygons[particle_num])-1,int(particle_arrays_pH[jj,particle_num])-1] += 1


    firstPldSwitch = False






    #print('num_good: {}'.format(num_good))
    #print('num_empty: {}'.format(num_empty))
    #print('num_bad: {}'.format(num_bad))

#    # Save the stationary settlers counts!
#    stationary_settlers_array = stationary_settlers_per_polygon
    
#    empty_trajectory_count_array = num_empty


#    end_script_time = time.time()
#    print(f'File processing took {(end_script_time-start_script_time)/60} minutes')
   
    # End of loop that used to be here
    # ------------------------------------------------------------------------------------------------
    #

#    # Remove stationary local settlers
#    num_stationary_settlers = 0
#    for polygon_num in range(num_polygons):
#        num_stationary_settlers_polygon = int(np.sum(stationary_settlers_array[polygon_num,:]))
#        #print(f'Number of stationary settlers: {num_stationary_settlers_polygon}')
#        prePdf_arrays_connectivity[array_dex,polygon_num,polygon_num] -= num_stationary_settlers_polygon
#        num_stationary_settlers += num_stationary_settlers_polygon
#        #prePdf_arrays_connectivity[array_dex,polygon_num,polygon_num] -= int(stationary_settlers_per_polygon[polygon_num])
#
#    print('Number of stationary settlers: {}'.format(num_stationary_settlers))









d = {}
d['particle_distances_per_pld_settlers'] = particle_distances_per_pld_settlers
d['particle_distances_per_pld_settlers_noStagnation'] = particle_distances_per_pld_settlers_noStagnation
d['particle_distances_per_pld_allTimesteps'] = particle_distances_per_pld_allTimesteps
d['particle_release_polygons'] = release_polygons
d['polygons_settled_per_particle'] = polygons_settled_per_particle


d["stagnant_settlers_per_particle_per_pld"] = stagnant_settlers_per_particle_per_pld
d["particle_settle_times_per_pld"] = particle_settle_times_per_pld
#particle_settle_times_per_pld_noStagnation



d['release_counts_per_polygon_array'] = release_counts_per_polygon_array
d['prePdf_arrays_connectivity'] = prePdf_arrays_connectivity
d['prePdf_arrays_connectivity_noStagnation'] = prePdf_arrays_connectivity_noStagnation
#    d['prePdf_arrays_settleTime'] = prePdf_arrays_settleTime
#    d['counter_array'] = counter_array
#    d['stationary_settlers_array'] = stationary_settlers_array
#    d['empty_trajectory_count_array'] = empty_trajectory_count_array
#d['pld_days'] = pld_list[pld_dex]
d['pld_list'] = pld_list
d['pld_list_timesteps'] = pld_list_timesteps

#    if time_dimension == 1:
#        d['prePdf_arrays_T'] = prePdf_arrays_T
    #d['prePdf_arrays_O2'] = prePdf_arrays_O2
    #d['prePdf_arrays_pH'] = prePdf_arrays_pH
    #d['O2_limit_list'] = O2_limit_list
    #d['pH_limit_list'] = pH_limit_list

save_output_file_name_pre = "conn_hist_data_"
save_output_file_name = save_output_file_name_pre + tracking_file.split('.')[0].split('/')[-1]
save_output_full_path = os.path.join(save_output_directory, save_output_file_name + f"_allPLDs_version_{script_version}_{polygon_string}")


np.savez(save_output_full_path, **d)

sample_start = 1000
sample_end = 1010

#print("polygons_settled_per_particle")
#print(polygons_settled_per_particle[sample_start:sample_end,:])
#print("particle_settle_times_per_pld")
#print(particle_settle_times_per_pld[sample_start:sample_end,:])
#print("stagnant_settlers_per_particle_per_pld")
#print(stagnant_settlers_per_particle_per_pld[sample_start:sample_end,:])
#print("particle_distances_per_pld_settlers")
#print(particle_distances_per_pld_settlers[sample_start:sample_end,:].astype(int))
###print("polygons_settled_per_particle")
##print(polygons_settled_per_particle[sample_start:sample_end,:])




