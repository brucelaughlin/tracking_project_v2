# Script to check if all opendrift jobs finished properly

import time
import numpy as np
import netCDF4
import argparse
import os
from pathlib import Path

t_start = time.time()

# Functions
def count_opendrift_jobs(filename,start_str,end_str):
    with open(filename,'r') as file:
        counting = False
        job_count = 0
        for line in file:
            if start_str in line:
                counting = True
            elif end_str in line:
                counting = False
            elif counting:
                job_count += 1
        return job_count

def get_config_params(filename):
#def get_config_params(filename,lifespan_str,numseed_str):
    with open(filename,'r') as file:
        for line in file:
            if lifespan_str in line:
                lifespan = int(line.split(':')[-1])
            elif numseed_str in line:
                numseed = int(line.split(':')[-1])
            elif dtsave_str in line:
                dtsave = int(line.split(':')[-1])
        return lifespan, numseed, dtsave

# Constant strings
config_str_1 = 'zstartNudgeList:'
config_str_2 = 'zznumberOfSeeds'
lifespan_str = 'drift:max_lifespan_days'
numseed_str = 'zznumberOfSeeds'
dtsave_str = 'runSave'
log_string_1 = 'active elements'
log_string_2 = 'number of floats seeded'
tracking_file_name_str = 'tracking_output_configFile_'
log_file_name_str = 'configFile_'
#report_directory = 'reports'
#report_name = 'run_report.txt'

parser = argparse.ArgumentParser()
parser.add_argument('rundir', type=str)
args = parser.parse_args()

run_dir = args.rundir

config_dir = os.path.join(run_dir, 'z_config_files')

config_dir_bytecode = os.fsencode(config_dir)

total_jobs_num = 0

for config_file in os.listdir(config_dir_bytecode):
    filename_pre = os.fsdecode(config_file)
    filename = os.path.join(config_dir,filename_pre)
    if filename.endswith(".yaml"):
        jobs_in_config = count_opendrift_jobs(filename, config_str_1, config_str_2)
        pld_max, num_seed_per_file, dt_save_minutes = get_config_params(filename)
        total_jobs_num += jobs_in_config

timesteps_output_per_day = 1440/dt_save_minutes


log_dir = os.path.join(run_dir, 'z_logs')



tracking_output_files = [filename for filename in os.listdir(run_dir) if filename.endswith(".nc")]
tracking_output_files.sort()

num_files = len(tracking_output_files)



#sample_dir = '/home/blaughli/test_dir'
#Path(os.path.join(sample_dir,report_directory).mkdir(exist_ok=True))        


print(run_dir)
print()
if num_files != total_jobs_num:
    print("Number of output files doesn't match the number of jobs specified in config files!")
else: 
    print('Correct number of output files, according to number jobs specified in config files')
print()

# I think the main thing is to compute the average trajectory length

sample_file_subname = 'configFile_013_job_13'

drift_timesteps_total = []

#sample_file_subname = 'configFile_013_job_13'

tracking_file_name_str = 'tracking_output_configFile_'
log_file_name_str = 'configFile_'

num_files = len(tracking_output_files)

for file_number in range(len(tracking_output_files)):
    
    tracking_output_file_pre = tracking_output_files[file_number]


    ###TESTING
    #if file_number > 3:
    #    break
    
    # Determine the number of floats in a single seeding
    #if file_number == 0:

    filename_split1 = tracking_output_file_pre.split(tracking_file_name_str)[-1]
    filename_split2 = filename_split1.split('_')
    config_file_num = filename_split2[0]
    filename_split3 = filename_split2[2]
    seed_num = filename_split3.split('.')[0]

    log_file = log_dir + '/' + log_file_name_str + config_file_num + '_job_' + seed_num + '.driftlog'
    
    with open(log_file,'r') as file:
        for line in file:
            if log_string_2 in line:
                floats_in_file = int(line.split(':')[-1])
                break
            else:
                floats_in_file = -9999

    
    tracking_output_file = os.path.join(run_dir,tracking_output_file_pre)
    dset = netCDF4.Dataset(tracking_output_file, 'r')
    lon_all = dset.variables['lon'][:]
    trajectory_timesteps_all = np.ma.MaskedArray.count(lon_all)
    average_trajectory_length = int(np.round((trajectory_timesteps_all/timesteps_output_per_day) / floats_in_file))


    num_max_life = 0
    num_min_life = 0
    timesteps_alive_list = []
    min_drift_time = pld_max

    for ii in range(np.shape(lon_all)[0]): 
        timesteps_alive_float = np.ma.MaskedArray.count(lon_all[ii,:])/timesteps_output_per_day
        timesteps_alive_list.append(timesteps_alive_float)
        drift_timesteps_total.append(timesteps_alive_float)
        if timesteps_alive_float == pld_max:
            num_max_life += 1
        if timesteps_alive_float == min_drift_time:
            num_min_life += 1
        if timesteps_alive_float < min_drift_time:
            min_drift_time = int(timesteps_alive_float)
            num_min_life = 1
    
    
    std_trajectory_length = int(np.round(np.std(np.array(timesteps_alive_list))))

    num_digits_max_count = len(str(floats_in_file))
    num_digits_max_drift = len(str(pld_max))
    num_digits_files = len(str(num_files))

    print(f"file {file_number+1:0{num_digits_files}}/{num_files}   Total floats: {floats_in_file}   Unit of time is days ->   PLD max: {pld_max},   Minimum drift time: {min_drift_time:0{num_digits_max_drift}},   n_max: {num_max_life:0{num_digits_max_count}},   n_min: {num_min_life:0{num_digits_max_count}},   average: {average_trajectory_length:03},   std: {std_trajectory_length:03}")



avg_trajectory_length_total = int(np.round(np.mean(np.array(drift_timesteps_total))))
std_trajectory_length_total = int(np.round(np.std(np.array(drift_timesteps_total))))

t_end = time.time()
report_runtime = t_end - t_start
report_runtime_minutes = int(np.round(report_runtime/60))
#report_runtime = int(np.round((t_start-t_end)/60))


print()
print()
print(f"Final statistics:   PLD max: {pld_max},   Average drift time: {avg_trajectory_length_total:03},   std drift time: {std_trajectory_length_total:03}")
print()
print()
print()
print(f"Report runtime: {np.round(report_runtime)} seconds = {report_runtime_minutes} minutes")
print(f"")
        

    
    
    




