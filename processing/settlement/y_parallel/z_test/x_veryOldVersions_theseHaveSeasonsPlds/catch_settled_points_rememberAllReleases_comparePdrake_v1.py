# rememberAllReleases - A new normalization to match Patrick's papers - we calulate the fraction of all released from a location that settled at a given
# location, rather then dividing only by the number of total settled

# v1: copied from v19 of the old normalization method



# create seasonal pdfs (djf, etc)

# Note that "status" is 0 when the particle is active, and a large magnitude negative
# number when not.  (strange!)


#---------------------------------------------------------------------
#---------------------------------------------------------------------
# I/O
#---------------------------------------------------------------------
#tracking_output_dir_pre = 'drift_150_physics_only_AKs_1en5_v1/'
#---------------------------------------------------------------------

# CHANGE THIS ONCE YOU'RE RE-ARRANGED THE OUTPUT DIRECTORY STRUCTURE
#tracking_output_base = "/home/blaughli/tracking_project/zz_output/tracking_project_output_projections/y_projections_1990_2020_KEEP/"
#tracking_output_base = "/data/blaughli/copiedFiles/y_projections_1990_2020_KEEP/"
#---------------------------------------------------------------------
#---------------------------------------------------------------------




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


#---------------------------------------------------------------------
#---------------------------------------------------------------------
# PARAMETERS and CONSTANTS

#parser = argparse.ArgumentParser()
#parser.add_argument("--trackingdir", type=str)
#parser.add_argument("--baseyear", type=int)
#args = parser.parse_args()
#tracking_output_dir = args.trackingdir
#base_year = args.baseyear

base_year = 1990
tracking_output_dir = '/data/blaughli/tracking_output/baseYear_1990/WC15N_GFDLTV_nRunsPerNode_15_nSeed_020_physicsOnly'


#---------------------------------------------------------------------
# Need to know the number of DAYS of each particle's life (fixed unless I change the 
# the deactivation time in the model code)
# NOTE: I believe that the extra point in time in all files is due to there being data saved at "time 0" - so,
# a '90 day run' with daily output will have 91 timesteps because initial state is saved.
#---------------------------------------------------------------------



# Create truncated normal distribution to draw from for settlement window openings
# Data taken from the report provided by Will (page 33, assuming mean is in the middle of PLD ranges, assuming ranges represent 3 stardard deviations)
#mu, sigma, lower, upper  = 120, 10, 90, 150


# Need to hardcode this for now, previous method of "calculating it" based on the input file doesn't work if we're
# using one input file for different maximum PLDs.
timesteps_per_day = 1


# WTF IS GOING ON HERE WITH THE TIMES???????

# Opendrift output times are seconds since Jan 1, 1979   (what??)
base_datetime = datetime.datetime(base_year,1,1,0,0,0)
#base_datetime = datetime.datetime(1988,1,1,0,0,0)


#---------------------------------------------------------------------
# Define the PLD for the species listed in Mallarie Yeager's report
#---------------------------------------------------------------------


## KELP BASS
pld_kelp_bass = [20,29]

## CALIFORNIA SHEEPHEAD
pld_ca_sheephead = [30,59]

## KELP ROCKFISH
pld_kelp_rockfish = [60,89]

## BLUE ROCKFISH, BLACK ROCKFISH
pld_blue_black_rockfish = [90, 149]


#pld_array=np.array([pld_kelp_bass,pld_ca_sheephead,pld_kelp_rockfish,pld_blue_black_rockfish])

#pld_array=np.array([pld_kelp_rockfish])

pld_array = np.array([[5,6],[10,11],[15,17],[20,22],[30,33],[45,49],[60,65],[90,98],[120,131],[150,164],[180,197]])

# choose one pld to compare, for now
pld_chosen_dex = 5

#---------------------------------------------------------------------
#---------------------------------------------------------------------


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
#---------------------------------------------------------------------

base_path = '/home/blaughli/tracking_project/'
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


tracking_output_files = [f for f in listdir(tracking_output_dir) if isfile(join(tracking_output_dir,f))]
tracking_output_files.sort()

# I always do this - is it bad practice?
tracking_output_dir = tracking_output_dir + "/"

save_output_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/'

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Prepare the data structure (2D array) for saving the pdf data

bounding_boxes_file_in = bounding_boxes_continent_dir + 'bounding_boxes_lonlat_coords_psi_coastline_wc15n_continent.p'
#bounding_boxes_file_in = bounding_boxes_continent_dir + 'bounding_boxes_lonlat_coords_{}_coastline_wc15n_continent.p'.format(point_type_line)
file = open(bounding_boxes_file_in,'rb')
boxes_lonlat = pickle.load(file)
file.close

n_boxes = 0
for box_lonlat in boxes_lonlat:
    n_boxes += 1

num_islands = 8
num_last_blob_island = 4
for island_dex in range(num_last_blob_island,num_islands+1):

        bounding_boxes_file_in = bounding_boxes_islands_dir + 'bounding_boxes_lonlat_wc15n_island_number_{}.p'.format(island_dex)

        file = open(bounding_boxes_file_in,'rb')
        boxes_lonlat = pickle.load(file)
        file.close        

        for box_lonlat in boxes_lonlat:
            n_boxes += 1
        #n_boxes += np.max([len(boxes_lonlat),len(boxes_lonlat[0])])



#---------------------------------------------------------------------
# This step is for setting up the pdf data structure which will be modified with each Opendrift output file.
# Also set up constants used throughout the runtime.

# number of pdfs per feild (4 seasonal and 1 overall = 5)
#n_pdfs = 5
n_pdfs = 1


# Need to have a way to see if I'm doing this right... without double counting errors, etc
tracking_output_file = tracking_output_dir + tracking_output_files[0]
dset = netCDF4.Dataset(tracking_output_file, 'r')
#status_all = dset.variables['status'][:]
particle_labels = dset.variables['trajectory'][:]
dset.close()
num_particles = len(particle_labels)
num_files = len(tracking_output_files)
counter_array=np.zeros((num_particles,num_files))

# For the histogram of average temperature experienced - just estimating range, since I don't know it without processing
#---------------------
# Make sure these match (0.1 = 1, 0.01 = 2, etc)
T_step = 0.1
n_decimals_round = 1
T_scale_factor = int(1/T_step)
#---------------------
T_min = 0
T_max = 30
n_T_steps = len(np.arange(T_min,T_max+1,T_step))



#---------------------------------------------------------------------
#---------------------------------------------------------------------
# hack for testing against patrick
# 1990 = nudges 3320 through 36400

tracking_output_files = []
for ii in range(3320,3680,40):
    tracking_output_files.append(tracking_output_dir + '/tracking_output_calcDT_060_saveDT_1440_buffer_100_nSeed_020_startNudge_{:06d}.nc'.format(ii))
num_files = len(tracking_output_files)
#---------------------------------------------------------------------
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# START THE MAIN LOOP!!!
#---------------------------------------------------------------------


for pld_dex in range(len(pld_array)):
        
    ###TESTING
    # test the 45-49 day pld
    if pld_dex != pld_chosen_dex:
        continue
        #break

    first_settlement_day = pld_array[pld_dex,0]
    last_settlement_day = pld_array[pld_dex,1]


    # THIS is the v5 adjustment: need to only use data from the pld, which begins after "first_settlement_day" and ends after "last_settlement_day"
    pld_length_days = last_settlement_day - first_settlement_day + 1


    timesteps_settlement_window = pld_length_days * timesteps_per_day
    timesteps_full_run = (last_settlement_day+1) * timesteps_per_day + 1


    first_settle_dex = first_settlement_day * timesteps_per_day +1
    last_settle_dex = (last_settlement_day+1) * timesteps_per_day + 1


    #---------------------------------------------------------------------
    # Create lists to store statistics for full run and seasonal subsets
    #---------------------------------------------------------------------

    # store the number of particles released from each box
    release_counts_per_cell = np.zeros((n_pdfs,n_boxes))


    # Connectivity (release box number vs settlement box number)
    pdf_arrays_connectivity = np.zeros((n_pdfs,n_boxes,n_boxes))

    # Time after PLD until settlement (saving only release location) (release box number vs settlement time)
    pdf_arrays_settleTime = np.zeros((n_pdfs,n_boxes,timesteps_settlement_window))

    # Number of days eposed to DO levels below 2.2 (saving only settlement location) (release box number vs exposure time)
    # note: Think the time dimension needs to be one bigger than the number of possible timesteps (ie need to include an option for "zero") 
    pdf_arrays_O2 = np.zeros((len(O2_limit_list),n_pdfs,n_boxes,timesteps_full_run))

    # Same idea for pH
    pdf_arrays_pH = np.zeros((len(pH_limit_list),n_pdfs,n_boxes,timesteps_full_run))

    # Histogram of average temperature experienced - just estimating range, since I don't know it without processing
    pdf_arrays_T = np.zeros((n_pdfs,n_boxes,n_T_steps))

    #---------------------------------------------------------------------
    #---------------------------------------------------------------------

    file_number = 0

    for tracking_output_file in tracking_output_files:
    #for tracking_output_file_pre in tracking_output_files:
       
        ###TESTING
        #if file_number > 9:
        #    break
        

            
        #tracking_output_file = tracking_output_dir + tracking_output_file_pre

        dset = netCDF4.Dataset(tracking_output_file, 'r')

        particle_labels = dset.variables['trajectory'][:]
        lon_all = dset.variables['lon'][:]
        lat_all = dset.variables['lat'][:]
        #z_all = dset.variables['z'][:]
        status_all = dset.variables['status'][:]
        time = np.array(dset['time'])
        # Exposure variables 
        O2_all = dset.variables['oxygen'][:]
        pH_all = dset.variables['pH'][:]
        temp_all = dset.variables['sea_water_temperature'][:]

        dset.close()

        O2_all *= conversion_factor


        # Prepare the list of possible seed months for the run
        run_seed_months_list = []
        for t in time:
            run_seed_months_list.append(datetime.datetime.strptime(str(base_datetime+datetime.timedelta(seconds=t)), '%Y-%m-%d %H:%M:%S').month)


        # Store the total number of particles
        num_particles = len(particle_labels)
          
        # Determine the output frequency, n per day.
        # Determine the timestep at which the settlement window opens
        trajectory_status = status_all[0,:]
        trajectory_mask = trajectory_status == 0




        # Create array to store all trajectory locations of all particles during their settlement window
        drift_lons = np.zeros((num_particles,timesteps_settlement_window))
        drift_lats = np.zeros((num_particles,timesteps_settlement_window))

        drift_O2 = np.zeros((num_particles,timesteps_full_run))
        drift_pH = np.zeros((num_particles,timesteps_full_run))
        drift_T = np.zeros((num_particles,timesteps_full_run))
        #drift_O2 = np.zeros((num_particles,timesteps_settlement_window))
        #drift_pH = np.zeros((num_particles,timesteps_settlement_window))
        #drift_T = np.zeros((num_particles,timesteps_settlement_window))

        starting_lons = []
        starting_lats = []
        seed_months = []

        #for particle_label in particle_labels:
        for particle_id in range(len(particle_labels)):
            #---------------------------------------------------------------------
            # indent here
            
            trajectory_status = status_all[particle_id,:]

            # Store month of first timestep of trajectory (ie seeding/starting month)
            seed_months.append(run_seed_months_list[np.where(trajectory_status == 0)[0][0]])

            trajectory_mask = trajectory_status == 0

               
            # Trim the trajectories to prepare for processing!
            particle_lon = lon_all[particle_id,trajectory_mask]
            particle_lat = lat_all[particle_id,trajectory_mask]
            particle_O2 = O2_all[particle_id,trajectory_mask]
            particle_pH = pH_all[particle_id,trajectory_mask]
            particle_T = temp_all[particle_id,trajectory_mask]

            # Get the STARTING coordinates of each particle
            starting_lons.append(particle_lon[0]) 
            starting_lats.append(particle_lat[0]) 
            
            ## Now that initial positions are known, cut out the "drifting window" and add what remains as a row
            ## in the <drift_lons>, <drift_lats> arrays
            #drift_lons[particle_id,:] = particle_lon[settlement_window_start:]
            #drift_lats[particle_id,:] = particle_lat[settlement_window_start:]

            #drift_O2[particle_id,:] = particle_O2[settlement_window_start:]
            #drift_T[particle_id,:] = particle_T[settlement_window_start:]
            
            drift_lons[particle_id,:] = particle_lon[first_settle_dex:last_settle_dex]
            drift_lats[particle_id,:] = particle_lat[first_settle_dex:last_settle_dex]
            drift_O2[particle_id,:] = particle_O2[0:last_settle_dex]
            drift_pH[particle_id,:] = particle_pH[0:last_settle_dex]
            drift_T[particle_id,:] = particle_T[0:last_settle_dex]
            #drift_O2[particle_id,:] = particle_O2[1:last_settle_dex]
            #drift_pH[particle_id,:] = particle_pH[1:last_settle_dex]
            #drift_T[particle_id,:] = particle_T[1:last_settle_dex]
            
            #drift_O2[particle_id,:] = particle_O2[settlement_window_start:settlement_window_end]
            #drift_pH[particle_id,:] = particle_pH[settlement_window_start:settlement_window_end]
            #drift_T[particle_id,:] = particle_T[settlement_window_start:settlement_window_end]

        # prepare lists to hold the starting/ending box numbers for each particle
        initial_boxes = np.zeros(num_particles).astype(int)
        settlement_boxes = np.zeros(num_particles).astype(int)
        # add settlement time storage
        settlement_times = np.zeros(num_particles).astype(int)

        # Exposure statistics
        particle_arrays_O2 = np.zeros((num_particles,len(O2_limit_list))).astype(int)
        particle_arrays_pH = np.zeros((num_particles,len(pH_limit_list))).astype(int)
        #particle_arrays_O2 = np.zeros((len(O2_limit_list),num_particles)).astype(int)
        #particle_arrays_pH = np.zeros((len(pH_limit_list),num_particles)).astype(int)  
        particle_array_T = np.zeros(num_particles)
        particle_array_driftTime = np.zeros(num_particles)


        #---------------------------------------------------------------------
        # FIRST MUST DETERMINE STARTING LOCATIONS
        #---------------------------------------------------------------------

        # create empty list to store all of the (starting) lat/lon pairs as tuples
        #points_lon_lat = []
        #for ii in range(num_particles):
        #    points_lon_lat.append((starting_lons[ii],starting_lats[ii]))
        #points_lon_lat = np.array(points_lon_lat)

        points_lon_lat = np.zeros((num_particles,2))
        for ii in range(num_particles):
            points_lon_lat[ii,0] = starting_lons[ii]
            points_lon_lat[ii,1] = starting_lats[ii]


        # Need a running "box_dex" index to keep track of box numbers across the calculations
        # Maybe don't have a "zero" box, so that 0 can mean NULL....
        box_dex = 0
        #box_dex = 1

        # WAIT... ISLANDS FIST, SO THEY'RE AT THE BOTTOM CORNER OF THE PDF PLOT...?!
        num_islands = 8
        num_last_blob_island = 4

        for island_dex in range(num_islands,num_last_blob_island-1,-1):
        #for island_dex in range(num_last_blob_island,num_islands+1):

            bounding_boxes_file_in = bounding_boxes_islands_dir + 'bounding_boxes_lonlat_wc15n_island_number_{}.p'.format(island_dex)

            # Load the boxes
            file = open(bounding_boxes_file_in,'rb')
            boxes_lonlat = pickle.load(file)
            file.close        

            #for box_lonlat in boxes_lonlat:
            for box_lonlat in reversed(boxes_lonlat):
                box_dex += 1
                if box_lonlat is not None:
                    path = plt_path.Path(np.transpose(box_lonlat))  # Transpose still needed?
                    particles_inside_flags = path.contains_points(points_lon_lat).astype(int)
                    initial_boxes = initial_boxes + (box_dex * particles_inside_flags)
                    #print(np.sum(particles_inside_flags))    


        # CONTINENT NEXT!!
        bounding_boxes_file_in = bounding_boxes_continent_dir + 'bounding_boxes_lonlat_coords_psi_coastline_wc15n_continent.p'
        #bounding_boxes_file_in = bounding_boxes_continent_dir + 'bounding_boxes_lonlat_coords_{}_coastline_wc15n_continent.p'.format(point_type_line)
        file = open(bounding_boxes_file_in,'rb')
        boxes_lonlat = pickle.load(file)
        file.close

        for box_lonlat in boxes_lonlat:
            box_dex += 1
            if box_lonlat is not None:
                path = plt_path.Path(np.transpose(box_lonlat))  # Transpose still needed?
                particles_inside_flags = path.contains_points(points_lon_lat).astype(int)
                initial_boxes = initial_boxes + (box_dex * particles_inside_flags)
                #box_dex += 1
                #print(np.sum(particles_inside_flags))    





        # NOW DETERMINE SETTLEMENT LOCATIONS!

        # create the "safety mask" that I'll use to make sure I only store the first box entered during the settlement window
        settlement_safety_mask = np.ones(num_particles).astype(int)
        #settlement_safety_mask = np.ones(num_particles)

        ## Sample the settlement window distribution for random settlement window opening times!
        #bio_windows_pre1 = biological_window_truncNorm_distribution.rvs(num_particles)
        #bio_windows = np.round(bio_windows_pre1)
        ## I think this actually belongs out here - mask needs to become more open with time
        #timestep_bio_window_mask=np.full(np.shape(bio_windows)[0],False)

        
        for time_dex in range(timesteps_full_run):
        #for time_dex in range(timesteps_settlement_window):

            print("file {}/{}, timestep {}/{}".format(file_number+1, num_files, time_dex,timesteps_full_run-1))
            #print("file {}/{}, timestep {}/{}".format(file_number+1, num_files, time_dex+1,timesteps_full_run))
            #print("file {}/{}, timestep {}/{}".format(file_number+1, num_files, time_dex+1,timesteps_settlement_window))

            #------------------------------------------------------
            # INDENT HERE

            # JUST FOR TESTING
            #time_dex = 0  # TESTING
        
            current_O2 = drift_O2[:,time_dex]
            current_pH = drift_pH[:,time_dex]
            current_T = drift_T[:,time_dex]

            # Exposure 
            for ii in range(len(O2_limit_list)):
                particle_arrays_O2[:,ii] = particle_arrays_O2[:,ii] + ((current_O2 < O2_limit_list[ii]).astype(int) * settlement_safety_mask)
                #particle_arrays_O2[ii,:] = particle_arrays_O2[ii,:] + ((current_O2 < O2_limit_list[ii]).astype(int) * settlement_safety_mask)
            for ii in range(len(pH_limit_list)):
                particle_arrays_pH[:,ii] = particle_arrays_pH[:,ii] + ((current_pH < pH_limit_list[ii]).astype(int) * settlement_safety_mask)
                #particle_arrays_pH[ii,:] = particle_arrays_pH[ii,:] + ((current_pH < pH_limit_list[ii]).astype(int) * settlement_safety_mask)
            particle_array_T = particle_array_T + (current_T * settlement_safety_mask)
            particle_array_driftTime = particle_array_driftTime + settlement_safety_mask
            
            

            #------------------------------------------------------
            #------------------------------------------------------
            # New modification: we want to track exposure for ENTIRE lives of larvae.  Hence the time condition on settlement computation
            #------------------------------------------------------
            if time_dex >= first_settle_dex:
                print("settle time dex: {}".format(time_dex))
            #------------------------------------------------------
            #------------------------------------------------------
            #------------------------------------------------------
            

                settle_window_dex = time_dex - first_settle_dex

                current_lons = drift_lons[:,settle_window_dex] 
                current_lats = drift_lats[:,settle_window_dex] 
                
                #points_lon_lat = []
                #for ii in range(num_particles):
                #    points_lon_lat.append((current_lons[ii],current_lats[ii]))
                #points_lon_lat = np.array(points_lon_lat)

                points_lon_lat = np.zeros((num_particles,2))
                for ii in range(num_particles):
                    points_lon_lat[ii,0] = current_lons[ii]
                    points_lon_lat[ii,1] = current_lats[ii]
                    #points_lon_lat[ii,0] = starting_lons[ii]
                    #points_lon_lat[ii,1] = starting_lats[ii]

                # Need a running "box_dex" index to keep track of box numbers across the calculations
                # Maybe don't have a "zero" box, so that 0 can mean NULL....
                box_dex = 0
                #box_dex = 1

                # WAIT... ISLANDS FIST, SO THEY'RE AT THE BOTTOM CORNER OF THE PDF PLOT...?!
                num_islands = 8
                num_last_blob_island = 4

                for island_dex in range(num_islands,num_last_blob_island-1,-1):
                #for island_dex in range(num_last_blob_island,num_islands+1):

                    bounding_boxes_file_in = bounding_boxes_islands_dir + 'bounding_boxes_lonlat_wc15n_island_number_{}.p'.format(island_dex)

                    # Load the boxes
                    file = open(bounding_boxes_file_in,'rb')
                    boxes_lonlat = pickle.load(file)
                    file.close        

                    #for box_lonlat in boxes_lonlat:
                    for box_lonlat in reversed(boxes_lonlat):
                        box_dex += 1
                        if box_lonlat is not None:
                            path = plt_path.Path(np.transpose(box_lonlat))  # Transpose still needed?
                            particles_inside_flags = path.contains_points(points_lon_lat)
                            particles_inside_flags_int = particles_inside_flags.astype(int)

                            settlement_boxes = settlement_boxes + (box_dex * particles_inside_flags_int) * settlement_safety_mask
                            settlement_times = settlement_times + ((settle_window_dex+1) * particles_inside_flags_int) * settlement_safety_mask
                            
                            # Update the counter of settled particles.  For checking consistency
                            counter_array[:,file_number] += particles_inside_flags_int * settlement_safety_mask
                            
                        
                            # Update the safety mask, so that settlement prevents further modifications to the stored settlement location
                            safety_mask_modifyer = ~particles_inside_flags
                            settlement_safety_mask = settlement_safety_mask * safety_mask_modifyer.astype(int)
                            


                # CONTINENT NEXT!!
                bounding_boxes_file_in = bounding_boxes_continent_dir + 'bounding_boxes_lonlat_coords_psi_coastline_wc15n_continent.p'
                #bounding_boxes_file_in = bounding_boxes_continent_dir + 'bounding_boxes_lonlat_coords_{}_coastline_wc15n_continent.p'.format(point_type_line)
                file = open(bounding_boxes_file_in,'rb')
                boxes_lonlat = pickle.load(file)

                for box_lonlat in boxes_lonlat:
                    box_dex += 1
                    if box_lonlat is not None:
                        path = plt_path.Path(np.transpose(box_lonlat))  # Transpose still needed?
                        particles_inside_flags = path.contains_points(points_lon_lat)
                        particles_inside_flags_int = particles_inside_flags.astype(int)

                        settlement_boxes = settlement_boxes + (box_dex * particles_inside_flags_int) * settlement_safety_mask
                        settlement_times = settlement_times + ((settle_window_dex+1) * particles_inside_flags_int) * settlement_safety_mask
                
                        # Update the counter of settled particles.  For checking consistency
                        counter_array[:,file_number] += particles_inside_flags_int * settlement_safety_mask
                            
                        # Update the safety mask, so that settlement prevents further modifications to the stored settlement location
                        safety_mask_modifyer = ~particles_inside_flags
                        settlement_safety_mask = settlement_safety_mask * safety_mask_modifyer.astype(int)
                        


        # Compute average temperature along trajectories
        average_T = [int(round(a/b, n_decimals_round) * T_scale_factor) for a,b in zip(particle_array_T, particle_array_driftTime)]



        # modify the pdf data structure
        for ii in range(num_particles):


            # Patrick only has the annual array
            array_dex_list = [0]

            # Need to skip the bad trajectories, which should have 0's in their initial boxes array
            if initial_boxes[ii] == 0:
                print('bad particle index : {}'.format(ii))
                continue
            else:

            #release_counts_per_cell[0,int(initial_boxes[ii])-1] += 1

            # wait we have 0's in here, so my logic requires we skip in this case
            #if initial_boxes[ii] == 0 or settlement_boxes[ii] == 0:
            #    if initial_boxes[ii] == 0:
            #        if drift_lons[ii,0] == 0:
            #            num_empty += 1
            #        elif drift_lons[ii,0] < extrema_lon[0] or drift_lons[ii,0] > extrema_lon[1] or drift_lats[ii,0] < extrema_lat[0] or  drift_lats[ii,0] > extrema_lat[1]:
            #            num_outofbounds += 1
            #        #elif drift_lons[ii,-1] == 0:
            #        #    num_empty += 1
            #    num_bad += 1
            #    continue
            #else:

            #    num_good += 1

                release_counts_per_cell[0,int(initial_boxes[ii])-1] += 1

                for array_dex in array_dex_list:
                    pdf_arrays_settleTime[array_dex,int(initial_boxes[ii])-1,int(settlement_times[ii])-1] += 1
                    pdf_arrays_connectivity[array_dex,int(initial_boxes[ii])-1,int(settlement_boxes[ii])-1] += 1

        file_number += 1


    d = {}
    d['release_counts_per_cell'] = release_counts_per_cell
    d['pdf_arrays_T'] = pdf_arrays_T
    d['pdf_arrays_O2'] = pdf_arrays_O2
    d['pdf_arrays_pH'] = pdf_arrays_pH
    d['pdf_arrays_connectivity'] = pdf_arrays_connectivity
    d['pdf_arrays_settleTime'] = pdf_arrays_settleTime
    d['counter_array'] = counter_array
    d['O2_limit_list'] = O2_limit_list
    d['pH_limit_list'] = pH_limit_list

    save_output_file_name_pre = "binned_data_seasonal_allReleases_baseYear_{}_".format(base_year)
    #save_output_file_name_pre = "pdf_data_seasonal_ranges_O2_pH_baseYear_{}_".format(base_year)
    save_output_file_name = save_output_file_name_pre + tracking_output_dir.split('/')[-2]
    save_output_full_path = save_output_directory + save_output_file_name + "_pld_{}_{}_forComparisonToPdrake".format(pld_array[pld_dex,0],pld_array[pld_dex,1])
    
    np.savez(save_output_full_path, **d)







