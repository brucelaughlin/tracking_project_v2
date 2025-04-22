# Script to check if all opendrift jobs finished properly

# Add on-the-fly functionality:
# - for unfinished runs, get % completed, and estimated time of/until completion


import sys
import datetime
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
#def get_config_params(filename,lifespan_str,numRunsJob_str):
    with open(filename,'r') as file:
        for line in file:
            if lifespan_str in line:
                lifespan = int(line.split(':')[-1])
            elif numRunsJob_str in line:
                numrunsjob = int(line.split(':')[-1])
            elif dtsave_str in line:
                dtsave = int(line.split(':')[-1])
            elif dtcalc_str in line:
                dtcalc = int(line.split(':')[-1])
        return lifespan, numrunsjob, dtsave, dtcalc

# Constant strings
config_str_1 = 'zstartNudgeList:'
config_str_2 = 'zznumberOfSeeds'
lifespan_str = 'drift:max_lifespan_days'
numRunsJob_str = 'numRunsPerJob:'
dtsave_str = 'runSave'
dtcalc_str = 'runCalc'
log_string_1 = 'active elements'
log_string_2 = 'number of floats seeded'
log_string_3 = '- step'
log_string_4 = 'INFO'
tracking_file_name_str = 'tracking_output_configFile_'
log_file_name_str = 'configFile_'
finished_str = 'output file size (raw)'
#report_directory = 'reports'
#report_name = 'run_report.txt'

minutes_in_day = 1440
seconds_in_day = 86400
seconds_in_hour = 3600

parser = argparse.ArgumentParser()
parser.add_argument('rundir', type=str)
parser.add_argument('numnodes', type=int)
args = parser.parse_args()

run_dir = args.rundir
num_nodes = args.numnodes


#config_dir_bytecode = os.fsencode(config_dir)

#config_file_list = os.listdir(config_dir_bytecode).sort()

config_dir = os.path.join(run_dir, 'z_config_files')
config_output_files = [filename for filename in os.listdir(config_dir) if filename.endswith(".yaml")]
config_output_files.sort()

total_jobs_num = 0
num_config = 0
num_runsPerJob_list = []

for config_file in config_output_files:
#for config_file in os.listdir(config_dir_bytecode):
    filename_pre = os.fsdecode(config_file)
    filename = os.path.join(config_dir,filename_pre)
    if filename.endswith(".yaml"):
        num_config += 1
        jobs_in_config = count_opendrift_jobs(filename, config_str_1, config_str_2)
        pld_max, num_runs_file, dt_save_minutes, dt_calc_minutes = get_config_params(filename)
        total_jobs_num += jobs_in_config
        num_runsPerJob_list.append(jobs_in_config)
        #num_runsPerJob_list.append(num_runs_file)

timesteps_output_per_save = minutes_in_day/dt_save_minutes


log_dir = os.path.join(run_dir, 'z_logs')


tracking_output_files = [filename for filename in os.listdir(run_dir) if filename.endswith(".nc")]
tracking_output_files.sort()

num_files = len(tracking_output_files)


# Printing parameters
num_digits_max_drift = len(str(pld_max))
num_digits_files = len(str(num_files))
num_digits_config = len(str(num_config))
num_digits_seed = len(str(int(max(num_runsPerJob_list))))
num_digits_floats = 0


drift_timesteps_total = []

tracking_file_name_str = 'tracking_output_configFile_'
log_file_name_str = 'configFile_'

config_summary_dict = {}

# Header information
print(run_dir)

any_finished_job_files = False
num_job_complete = 0
firstFile = True
floats_in_file = -99999999
report_string_list = []
usable_log_exists = False
config_dex = 0

for file_number in range(len(tracking_output_files)):
    
    tracking_output_file_pre = tracking_output_files[file_number]

    ###TESTING
    #if file_number > 3:
    #    break
    
    # Determine the number of floats in a single seeding
    #if file_number == 0:

    filename_split1 = tracking_output_file_pre.split(tracking_file_name_str)[-1]
    filename_print = tracking_output_file_pre
    #filename_print = filename_split1
    filename_split2 = filename_split1.split('_')
    config_file_num_str = filename_split2[0]
    if np.logical_not(firstFile):
        if config_file_num_str != config_file_num_str_prevoius:
            config_dex += 1
    config_file_num_str_prevoius = config_file_num_str
    filename_split3 = filename_split2[2]
    job_num_str = filename_split3.split('.')[0]

    log_file = log_dir + '/' + log_file_name_str + config_file_num_str + '_job_' + job_num_str + '.driftlog'
   
    timestep_current = 0
    timestep_final = 0
    num_timesteps = 0
    time_final = 0

    finished_switch = False
    first_timestep_switch = True

    num_full_days_run_floor = 0
    time_previous = 0

    time_per_complete_day_list = []

    with open(log_file,'r') as file:
        for line in file:
            if finished_str in line:
                finished_switch = True
                any_finished_job_files = True
                num_job_complete += 1
            if log_string_2 in line:
                if firstFile:
                    floats_in_file = int(line.split(':')[-1])
                    num_digits_floats = len(str(floats_in_file))
                    firstFile = False
            if log_string_3 in line:
                num_timesteps += 1
                timestep_info_pre1 = line.split(log_string_3)[-1]
                timestep_info_pre2 = timestep_info_pre1.split(' - ')[0]
                timestep_info = timestep_info_pre2.split('of')
                timestep_current = int(timestep_info[0])
                timestep_final = int(timestep_info[1])
                if first_timestep_switch:
                    timestep_info_pre1 = line.split(log_string_4)[-0]
                    timestep_info = timestep_info_pre1.split(':')
                    time_h = int(timestep_info[0])
                    time_m = int(timestep_info[1])
                    time_s = int(timestep_info[2])
                    time_initial = time_h *3600 + time_m * 60 + time_s
                    first_timestep_switch = False
                    time_previous = time_initial

                else:
                    usable_log_exists = True
                    timestep_info_pre1 = line.split(log_string_4)[-0]
                    timestep_info = timestep_info_pre1.split(':')
                    time_h = int(timestep_info[0])
                    time_m = int(timestep_info[1])
                    time_s = int(timestep_info[2])
                    time_final = time_h *3600 + time_m * 60 + time_s
                    # Account for crossing into next day
                    if time_final < time_previous:
                        total_time_previous_day = (24*60*60 - time_previous) + (time_previous - time_initial)
                        time_per_complete_day_list.append(total_time_previous_day)
                        time_initial = 0
                        #num_full_days_run_floor += 1


                    time_previous = time_final
        
    # Need to handle case where current file only has one timestep written - in this case, can't generate estimage, so skip the file
    if time_final != 0:
        runtime_hrs = (sum(time_per_complete_day_list) + time_final - time_initial)/seconds_in_hour
   

        if config_file_num_str not in config_summary_dict:
            config_summary_dict[config_file_num_str] = {"runtime":[],"fraction_complete":[],"finished":[]}

        config_summary_dict[config_file_num_str]['runtime'].append(runtime_hrs)
        config_summary_dict[config_file_num_str]['fraction_complete'].append(num_timesteps/timestep_final)
        config_summary_dict[config_file_num_str]['finished'].append(finished_switch)


        if finished_switch:

            tracking_output_file = os.path.join(run_dir,tracking_output_file_pre)
            dset = netCDF4.Dataset(tracking_output_file, 'r')
            lon_all = dset.variables['lon'][:]
            trajectory_timesteps_all = np.ma.MaskedArray.count(lon_all)
            average_trajectory_length = int(np.round((trajectory_timesteps_all/timesteps_output_per_save) / floats_in_file))

            num_max_life = 0
            num_min_life = 0
            timesteps_alive_list = []
            min_drift_time = pld_max

            for ii in range(np.shape(lon_all)[0]): 
                timesteps_alive_float = np.ma.MaskedArray.count(lon_all[ii,:])/timesteps_output_per_save
                timesteps_alive_list.append(timesteps_alive_float)
                drift_timesteps_total.append(timesteps_alive_float)
                # BIG ISSUE HERE - I NEED TO UNDERSTAND WHY I'M GETTING MORE OUTPUT TIMESTEPS THAN THE MAX PLD WOULD INDICATE (I GET 2 EXTRA DAYS, NOT SURE IF THAT'S ALWAYS THE CASE)
                if timesteps_alive_float >= pld_max:
                #if timesteps_alive_float == pld_max:
                    num_max_life += 1
                if timesteps_alive_float == min_drift_time:
                    num_min_life += 1
                if timesteps_alive_float < min_drift_time:
                    min_drift_time = int(timesteps_alive_float)
                    num_min_life = 1
            
            
            std_trajectory_length = int(np.round(np.std(np.array(timesteps_alive_list))))

            report_string = f"configFile {config_dex+1}/{num_config},   {filename_print},    n_max_drift: {num_max_life:0{num_digits_floats}},   Minimum drift time: {min_drift_time:0{num_digits_max_drift}},   n_min_drift: {num_min_life:0{num_digits_floats}},   average: {average_trajectory_length:03},   std: {std_trajectory_length:03},   percent of max timesteps used: {100*num_timesteps/timestep_final:02.2f}"
            
            #report_string = f"configFile {config_dex+1}/{num_config},   configFile_job: {filename_print},    n_max_drift: {num_max_life:0{num_digits_floats}},   Minimum drift time: {min_drift_time:0{num_digits_max_drift}},   n_min_drift: {num_min_life:0{num_digits_floats}},   average: {average_trajectory_length:03},   std: {std_trajectory_length:03},   percent of max timesteps used: {100*num_timesteps/timestep_final:02.2f}"
            
            #report_string = f"config file {int(config_file_num_str):0{num_digits_config}}/{num_config-1:0{num_digits_config}}, config job {int(job_num_str):0{num_digits_seed}}/{num_runsPerJob_list[config_dex]-1:0{num_digits_seed}}  ->  n_max_drift: {num_max_life:0{num_digits_floats}},   Minimum drift time: {min_drift_time:0{num_digits_max_drift}},   n_min_drift: {num_min_life:0{num_digits_floats}},   average: {average_trajectory_length:03},   std: {std_trajectory_length:03},   percent of max timesteps used: {100*num_timesteps/timestep_final:02.2f}"
            
            report_string_list.append(report_string)

# Exit with error message if every log file only contains one timestep record
if np.logical_not(usable_log_exists):
    sys.exit("Too early to generate report; Need at least two completed timesteps")



#print()
#print()
#if num_files != total_jobs_num:
#    print("Number of output files doesn't match the number of jobs specified in config files - run is incomplete, or an error/mistake has occurred")
#else: 
#    print('Correct number of output files, according to number jobs specified in config files')


#config_summary_dict[config_file_num_str]['runtime'].append(runtime_hrs)
#config_summary_dict[config_file_num_str]['fraction_complete'].append(num_timesteps/timestep_final)
#config_summary_dict[config_file_num_str]['finished'].append(finished_switch)


# Predict overall runtime, and time of finish:
any_finished_config_files = False
num_config_complete = 0
runtime_avg_pre = 0
num_partial = 0
fraction_avg_pre = 0
runtime_partial_avg_pre = 0
partial_time_left_avg_pre = 0

for config_file_num in config_summary_dict.keys():
    if sum(config_summary_dict[config_file_num]['finished']) == len(config_summary_dict[config_file_num]['finished']):
        any_finished_config_files = True
        num_config_complete += 1
        runtime_avg_pre += sum(config_summary_dict[config_file_num]['runtime'])/len(config_summary_dict[config_file_num]['runtime'])
    else:
        num_partial += 1
        fraction_avg_internal = sum(config_summary_dict[config_file_num]['fraction_complete'])/len(config_summary_dict[config_file_num]['fraction_complete'])
        #fraction_avg_pre += fraction_avg_internal
        runtime_partial_avg_internal = sum(config_summary_dict[config_file_num]['runtime'])/len(config_summary_dict[config_file_num]['runtime'])
        partial_time_left_avg_pre += runtime_partial_avg_internal * (1/fraction_avg_internal - 1)
        #partial_time_left_avg_pre += runtime_partial_avg_internal/fraction_avg_internal


# Overall statistics:
if any_finished_config_files:
    avg_trajectory_length_total = int(np.round(np.mean(np.array(drift_timesteps_total))))
    std_trajectory_length_total = int(np.round(np.std(np.array(drift_timesteps_total))))
    runtime_avg = runtime_avg_pre/num_config_complete

    

#runtime_partial_avg = runtime_partial_avg_pre/num_partial
#fraction_avg = fraction_avg_pre/num_partial
if num_partial > 0:
    partial_time_left_avg = partial_time_left_avg_pre/num_partial 
else:
    partial_time_left_avg = 0

#print()
#print(runtime_avg_pre)
#print(f"Number of completed config files {num_config_complete}")
#print(f"Average completed file runtime:       {runtime_avg}")
#print(f"Average incomplete estimated runtime: {partial_time_left_avg}")



unprocessed_jobs_runtime = 0

num_config_unprocessed = num_config - num_config_complete - num_partial

partial_jobs_runtime = partial_time_left_avg * num_partial 

if any_finished_config_files:
    unprocessed_jobs_runtime = runtime_avg * (np.ceil(num_config_unprocessed/num_nodes))
else:
    unprocessed_jobs_runtime = partial_time_left_avg * (np.ceil(num_config_unprocessed/num_nodes))

time_required_hours = unprocessed_jobs_runtime + partial_jobs_runtime

time_now = datetime.datetime.now()
time_estimated_completion = time_now + datetime.timedelta(hours=time_required_hours)

ttc_pre = datetime.timedelta(hours=time_required_hours)


ttc = {"days": ttc_pre.days}
ttc["hours"], rem = divmod(ttc_pre.seconds, 3600)
ttc["minutes"], ttc["seconds"] = divmod(rem, 60)


# String formatting stuff
#if ttc["days"] > 0:
#    day_str = "days"
#    if ttc["days"] > 1:
#        day_str = "day"
#    if ttc["minutes"] > 30:
#        ttc["hours"] += 1
#if ttc["hours"] > 0:
#    hour_str = "hours"
#    if ttc["days"] > 0:
#else:
#    hour_str = "hour"
#if ttc["minutes"] > 0:
#    minute_str = "minutes"
#else:
#    minute_str = "minute"
#if ttc["seconds"] > 0:
#    second_str = "seconds"
#else:
#    second_str = "second"
#

if ttc["days"] > 0:
    fmt = "days: {days}, hours: {hours}"
elif ttc["hours"] > 0:
    fmt = "hours: {hours}, minutes: {minutes}"
else:
    fmt = "minutes: {minutes}, seconds: {seconds}"

#if ttc["days"] > 0:
#    if ttc["days"] > 1:
#        fmt = "{days} days {hours} hours"
#    else:
#        fmt = "{days} day {hours} hours"
#elif ttc["hours"] > 1:
#    fmt = "{hours} hours {minutes} minutes"
#else:
#    fmt = "{seconds} hours {seconds} seconds"

ttc_string = fmt.format(**ttc)

time_now_print = time_now.strftime('%m-%d %H:%M')
time_estimated_completion_print = time_estimated_completion.strftime('%m-%d %H:%M')

if num_config_complete < num_config:
    print()
    print(f"Config files: {num_config_unprocessed} unprocessed, {int(num_partial)} processing, {num_config_complete} completed")
    #print(f"Config file processing status: {num_config_unprocessed} unprocessed, {int(num_partial)} currently being processed, {num_config_complete} completed")

    print()
    #print()
    #print(f'Current time:       {time_now_print}')
    #print()
    #print()
    print(f'Time remaining:     {ttc_string}')
    print(f'Completion time~:   {time_estimated_completion_print}')
    if np.logical_not(any_finished_config_files):
        print()
        print("Note: Zero config files have run to completion; time estimates will improve after all jobs from at least one config file have finished.")

    


else:
    
    trt_pre = datetime.timedelta(hours=num_config * runtime_avg)
    trt = {"days": trt_pre.days}
    trt["hours"], rem = divmod(trt_pre.seconds, 3600)
    trt["minutes"], trt["seconds"] = divmod(rem, 60)
    if trt["minutes"] > 30:
        trt["hours"] += 1

    if trt["days"] > 0:
        fmt = "days: {days}, hours: {hours}"
    elif trt["hours"] > 0:
        fmt = "hours: {hours}, minutes: {minutes}"
    else:
        fmt = "minutes: {minutes}, seconds: {seconds}"

    trt_string = fmt.format(**trt)
    print()
    print(f'Run complete, total runtime: {trt_string}')
    #print('Run complete!')


if any_finished_config_files:
    print()
    print(f"Overall statistics for completed runs:   PLD max: {pld_max},   Average drift time: {avg_trajectory_length_total:03},   std drift time: {std_trajectory_length_total:03}")


    #print()
    print()
    print('Data for completed runs:')
    print()
    for report_string in report_string_list:
        print(report_string)




#t_end = time.time()
#report_runtime = t_end - t_start
#report_runtime_minutes = int(np.round(report_runtime/60))
#
#print()
#print()
#
#if report_runtime_minutes > 1:
#    print(f"Report runtime: {report_runtime_minutes} minutes")
#else:
#    print(f"Report runtime: {report_runtime_minutes} minute")
        

    
    
    




