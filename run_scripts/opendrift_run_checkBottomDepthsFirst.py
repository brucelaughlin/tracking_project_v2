
import yaml
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

import logging
import subprocess
import pickle
import netCDF4
import numpy as np
import datetime
from datetime import timedelta
from dateutil.relativedelta import relativedelta
import time
import sys 
import os
from pathlib import Path
import argparse
#from opendrift_custom.models.larvaldispersal_track_eco_variables_RecordKicksTest import LarvalDispersal
from opendrift_custom.models.larvaldispersal_track_eco_variables import LarvalDispersal

logger = logging.getLogger('opendrift_run_v2')



# Track how long this takes to run
t_init = time.time()


parser = argparse.ArgumentParser()
parser.add_argument("--configfile", type=str, help='help yourself')
parser.add_argument("--level", default='DEBUG', type=str, help='debug level string, ie INFO, WARNING, DEBUG,..')
parser.add_argument("--jobrunnumber", type=int, help='blah')
parser.add_argument("--jobrunstring", type=str, help='blah')

args = parser.parse_args()

# paul's tricky one line in-line dictionary hack
try:
    logging.basicConfig(level={'info':logging.INFO,'debug':logging.DEBUG,'warning':logging.WARNING,'error':logging.ERROR}[args.level.lower()])
except KeyError:
    parser.error('Invalid --level string')

config_file = args.configfile
job_run_number = args.jobrunnumber
job_run_string = args.jobrunstring

stream = open(config_file,'r')
#config_dict = yaml.load(stream, Loader=yaml.BaseLoader)    
config_dict = yaml.safe_load(stream)    
stream.close()

#polygon_seed_lon_lat_file = config_dict['polygonSeedLonLatFile']
#polygon_seed_ii_jj_file = config_dict['polygonSeedIJFile']
polygon_seed_lon_lat_file = os.path.join(config_dict['seedLocationDataDir'],"seed_coordinates_lon_lat.p")
polygon_seed_ii_jj_file = os.path.join(config_dict['seedLocationDataDir'],"seed_coordinates_ii_jj.p")

roms_forcing_switch = config_dict["romsForcingSwitch"]
if roms_forcing_switch == 'true':
    roms_forcing_switch = True
elif roms_forcing_switch == 'false':
    roms_forcing_switch = False

if roms_forcing_switch:
    from opendrift.readers.reader_ROMS_native import Reader
else:
    from opendrift.readers.reader_netCDF_CF_generic import Reader


grid_path = config_dict["gridPath"]
run_calc = config_dict["runCalc"]
run_save = config_dict["runSave"]
buffer_length_export = config_dict["bufferLengthExport"]
number_of_seeds = config_dict["zznumberOfSeeds"]
days_between_seeds = config_dict["seedSpacing"]
base_year = config_dict["baseYear"]
his_dir = config_dict["inputDir"]
output_dir = config_dict["outputDir"]
start_nudge = config_dict["zstartNudgeList"][job_run_number]


run_details_string = 'calcDT_{b:03d}_saveDT_{c:04d}_exportBuffer_{d:03d}_nSeed_{e:03d}_startNudge_{f:06d}'.format(b=run_calc,c=run_save,d=buffer_length_export,e=number_of_seeds,f=start_nudge)

print('USER PRINT STATEMENT: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',flush=True)
print('USER PRINT STATEMENT: {}'.format(run_details_string),flush=True)
print('USER PRINT STATEMENT: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',flush=True)
#-------------------------------------------------


# -----------------------------------------------------------------------------
# run configuration parameters:
seed_window_length = (number_of_seeds - 1) * days_between_seeds + 1
save_dt = run_save * 60;
run_dt = run_calc * 60

base_datetime = datetime.datetime(base_year,1,1,12,0,0)

# This is dumb
his_file_wildcard = "*.nc"
#his_file_wildcard = '*.nc'
#his_file_wildcard = 'wc15*.nc'
# -----------------------------------------------------------------------------

try:
    if os.path.splitext(grid_file)[-1] == ".npz":
        d = np.load(grid_file)
        h = d['h']
    elif os.path.splitext(grid_file)[-1] == ".nc":
        d = netCDF4.Dataset(grid_file)
        h = np.array(d['h'])
        d.close()
    else:
        raise ValueError("Grid files are expected to have either '.nc' or '.npz' as extensions")
except ValueError as e:
    print(f"Error: {e}")


#-------- History Files -----------------
his_files = os.path.join(his_dir,his_file_wildcard)

#----------Output netCDF File---------------------
tracking_output_pre = 'tracking_output_{}.nc'.format(job_run_string)
tracking_output_file = os.path.join(output_dir,tracking_output_pre)

# prepare memory plot file name
output_file_split = tracking_output_file.split('.')
output_file_pre = output_file_split[0]
output_png_file = output_file_pre + '.png'

#--------- Get seeding lon/lat coordinates data ------------
file = open(polygon_seed_lon_lat_file,'rb')
list_of_arrays_of_points_in_polygons_lonlat = pickle.load(file)
file.close
file = open(polygon_seed_ii_jj_file,'rb')
list_of_arrays_of_points_in_polygons_iijj = pickle.load(file)
file.close



# We want profiles of floats at each lat/lon starting point.  Space them "depth_step" (5m) apart, from the surface
# down to a set depth min "min_float_depth" (20m) or the bottom depth, whichever is shallower
min_float_depth = 20
#min_float_depth = 15
depth_step = 5


# Compute the seeding time window
start_seed_time = base_datetime + relativedelta(days = start_nudge)
end_seed_time = start_seed_time + relativedelta(days = seed_window_length)


# Prepare lists of coordinates for seeding
lons = []
lats = []
zs = []
times = []

test_switch_horizontal = config_dict["testSwitchHorizontal"]
if test_switch_horizontal == 'true':
    test_switch_horizontal = True
elif test_switch_horizontal == 'false':
    test_switch_horizontal = False

test_switch_vertical = config_dict["testSwitchVertical"]
if test_switch_vertical == 'true':
    test_switch_vertical = True
elif test_switch_vertical == 'false':
    test_switch_vertical = False

# Note that I am changing over to using only CSV files as input for the polygon coordinates, and 
# this changes the order of the dimensions, I believe.  So, for now I'm commenting out the old
# code but leaving it in case I mess something up


if test_switch_horizontal:
    test_cells = [179,180,181,182] # Mercator, coastal cell
    #test_cells = [181] # Mercator, coastal cell
    #test_cells = [32] # WC15
    #test_cells = [26,27,29,32]
    for run_day in range(0,seed_window_length,days_between_seeds):
        for ii in test_cells:
            for jj in range(np.shape(list_of_arrays_of_points_in_polygons_lonlat[ii])[1]):
                depth_min = min_float_depth
                for kk in range(1):
                    zs.append(-kk*depth_step)
                    lons.append(list_of_arrays_of_points_in_polygons_lonlat[ii][0,jj])
                    lats.append(list_of_arrays_of_points_in_polygons_lonlat[ii][1,jj])
                    times.append(datetime.datetime.strptime(str(start_seed_time+datetime.timedelta(days=run_day)), '%Y-%m-%d %H:%M:%S'))






elif test_switch_vertical:
    #lat_list = [36.2, 36.2]
    #lon_list = [-121.176, -123]
    #lat_list = [35.5333] 
    #lon_list = [-121.867]
    lat_list = [36.0667] 
    lon_list = [-122.5]
    #lon_list = [-121.967]
    ## h at this location is 1050.34
        
    opendrift_config_dict = config_dict['modelConfigDict']
    target_depth = -1 * opendrift_config_dict['drift:target_depth_day']

    num_particles_per_xy_location = 1000

    #test_cell = 26
    run_day = 0
    ###for jj in range(np.shape(list_of_polygon_vertex_lonlat_arrays[test_cell])[1]):
    #for jj in range(1):
     
        #bottom_depth = h[points_in_boxes_i_j[test_cell][0,jj],points_in_boxes_i_j[test_cell][1,jj]]
        #depth_min = np.floor(min(min_float_depth,bottom_depth))

    ### offshore particle
    #for ii in range(1,2)):
    
    ### nearshore particle
    #for ii in range(1)):

    ### both particles
    for ii in range(len(lon_list)):
        for jj in range(num_particles_per_xy_location):
            #zs.append(target_depth + np.random.uniform(target_depth - 500, target_depth + 500))
            zs.append(target_depth)
            #zs.append(-kk*depth_step)
            lons.append(lon_list[ii])
            lats.append(lat_list[ii])
            times.append(datetime.datetime.strptime(str(start_seed_time+datetime.timedelta(days=run_day)), '%Y-%m-%d %H:%M:%S'))
   
#    # seed the same number of particles at half of the target depth (to get motion unaffected by bottom, ideally)
#    for ii in range(len(lon_list)):
#        for jj in range(num_particles_per_xy_location):
#            zs.append(target_depth*0.5)
#            #zs.append(-kk*depth_step)
#            lons.append(lon_list[ii])
#            lats.append(lat_list[ii])
#            times.append(datetime.datetime.strptime(str(start_seed_time+datetime.timedelta(days=run_day)), '%Y-%m-%d %H:%M:%S'))


else:
    for run_day in range(0,seed_window_length,days_between_seeds):
        for ii in range(len(list_of_arrays_of_points_in_polygons_lonlat)):
            for jj in range(np.shape(list_of_arrays_of_points_in_polygons_lonlat[ii])[1]):
                bottom_depth = h[list_of_arrays_of_points_in_polygons_iijj[ii][0,jj],list_of_arrays_of_points_in_polygons_iijj[ii][1,jj]]
                depth_min = np.floor(min(min_float_depth,bottom_depth))
                for kk in range(int(np.floor(depth_min / depth_step)) + 1):
                    zs.append(-kk*depth_step)
                    lons.append(list_of_arrays_of_points_in_polygons_lonlat[ii][0,jj])
                    lats.append(list_of_arrays_of_points_in_polygons_lonlat[ii][1,jj])
                    times.append(datetime.datetime.strptime(str(start_seed_time+datetime.timedelta(days=run_day)), '%Y-%m-%d %H:%M:%S'))




print('USER PRINT STATEMENT: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',flush=True)
print('USER PRINT STATEMENT: number of floats seeded: {} '.format(len(lons)),flush=True)
print('USER PRINT STATEMENT: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',flush=True)

lons = np.asarray(lons)
lats = np.asarray(lats)
zs = np.asarray(zs)


# ------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------
# Initialize Opendrift
o = LarvalDispersal(loglevel=config_dict['logLevel'])

print('USER PRINT STATEMENT: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',flush=True)
# Execute configuration directives using our config file
opendrift_config_dict = config_dict['modelConfigDict']

### Turn off mixing if vertical swimming is activated
#if opendrift_config_dict['drift:vertical_swim_speed'] > 0:
#    opendrift_config_dict['drift:vertical_mixing'] = False 

for key,value in opendrift_config_dict.items():
    # yaml dumping is turning booleans within dictionaries into lowercase words!!!
    if value == 'true':
        value = True
    elif value == 'false':
        value = False

    o.set_config(key,value)
    print(f'USER PRINT STATEMENT: setting o.set_config({key}, {value})',flush=True)
print('USER PRINT STATEMENT: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',flush=True)

particle_lifetime = opendrift_config_dict['drift:max_lifespan_days']
# Calculate run duration in hours
duration_safety_buffer = 3
run_duration_hours = (seed_window_length + particle_lifetime + duration_safety_buffer) * 24

# ------------------------- Adding Reader --------------------------------
t_read_0 = time.time()

new_variables = config_dict['newVariables']

#reader_list = []

# It's possible that user specified no new variables in the config file
if new_variables is None:
    reader_current_job = Reader(his_files)
else:
    reader_current_job = Reader(his_files, standard_name_mapping=new_variables)
o.add_reader(reader_current_job)

t_read_1 = time.time()
reader_time = t_read_1 - t_read_0

print('USER PRINT STATEMENT: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',flush=True)
print('USER PRINT STATEMENT: Minutes taken to add readers for 3 years: {}'.format(round(reader_time)/3600,3),flush=True)
print('USER PRINT STATEMENT: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',flush=True)
#--------------------------------------------------------------------------

export_variables = config_dict['exportVariables']

t_run_start = time.time()

o.seed_elements(lon=lons,lat=lats, z=zs, time=times, origin_marker = 0)

o.run(duration=timedelta(hours=run_duration_hours), time_step=run_dt, time_step_output=save_dt, outfile = tracking_output_file, export_variables = export_variables, export_buffer_length=buffer_length_export)

# plot memory usage, safe figure
#o.plot_memory_usage(filename=output_png_file)

t_run_end = time.time()
total_runtime = t_run_end-t_run_start
total_execution_time = t_run_end-t_init

summary_string = '{}, time_start: {}, time_end: {},  number_of_seeds: {}, days_particle_lifetime: {}, readerloading (mins): {}, run_time (hrs): {}, execution_time (hrs): {}\n'.format(run_details_string,str(t_init),str(t_run_end),str(number_of_seeds),particle_lifetime,round(reader_time/60,3),round(total_runtime/3600,3), round(total_execution_time/3600,3))

print('USER PRINT STATEMENT: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',flush=True)
print('USER PRINT STATEMENT: ' + str(o),flush=True)
#print('USER PRINT STATEMENT: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',flush=True)
#print('USER PRINT STATEMENT: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',flush=True)
print('USER PRINT STATEMENT: \nprogram start time: {}\n'.format(t_init),flush=True)
print('USER PRINT STATEMENT: \nopendrift start time: {}\n'.format(t_init),flush=True)
print('USER PRINT STATEMENT: \nend time: {}\n'.format(t_run_end),flush=True)
#print('USER PRINT STATEMENT: \ntotal runtime: {} hours\n'.format(round(total_runtime/3600,1)),flush=True)
#print('USER PRINT STATEMENT: \ntotal execution time: {} hours\n'.format(round(total_execution_time/3600,1)),flush=True)
print('USER PRINT STATEMENT: \ntotal runtime: {} hours\n'.format(round(total_runtime/3600,4)),flush=True)
print('USER PRINT STATEMENT: \ntotal execution time: {} hours\n'.format(round(total_execution_time/3600,4)),flush=True)
print('USER PRINT STATEMENT: \nsummary info: {}\n'.format(summary_string),flush=True)
#print('USER PRINT STATEMENT: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',flush=True)


        

# Compress the output file
if np.logical_not(test_switch_horizontal or test_switch_vertical):
#if np.logical_not(test_switch_horizontal) or np.logical_not(test_switch_vertical):
    
    bash_command = "ls -lh {}".format(tracking_output_file)
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    size_raw = process.stdout.read()

    bash_command = "nc_compress {}".format(tracking_output_file)
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    time.sleep(120)

    bash_command = "ls -lh {}".format(tracking_output_file)
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    size_compressed = process.stdout.read()

    #print('USER PRINT STATEMENT: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv',flush=True)
    print('USER PRINT STATEMENT: \noutput file size (raw): {}\n'.format(size_raw),flush=True)
    print('USER PRINT STATEMENT: \noutput file size (compressed): {}\n'.format(size_compressed),flush=True)
    #print('USER PRINT STATEMENT: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',flush=True)




print('Finished',flush=True)
print('USER PRINT STATEMENT: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',flush=True)





