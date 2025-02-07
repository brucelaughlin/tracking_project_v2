# Generate config files for Opendrift Runs

# Now using the yearly files, not daily within each year

import yaml
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

import argparse
import shutil
import os
from pathlib import Path
import numpy as np
import subprocess
import re
import pdb

config_template_file = '/home/blaughli/tracking_project_v2/config_files_for_setup/z_templates/config_template.config.yaml'

tracking_dir_pre_str = 'tracking_run_'

parser = argparse.ArgumentParser()
parser.add_argument('configfile', type=str)
parser.add_argument('testSwitch', nargs='?', type=str)
args = parser.parse_args()

config_for_config_file = args.configfile

stream = open(config_for_config_file,'r')
#config_dict = yaml.load(stream, Loader=yaml.BaseLoader)
config_dict = yaml.safe_load(stream)
stream.close()

his_file_name_pre = config_dict['hisFileNamePre']
base_year = config_dict['baseYear']
inputDir = config_dict['inputDir']
baseOutputDir = config_dict['baseOutputDir']
outputDirName = config_dict['outputDirName']

numRunsPerJob = config_dict['numRunsPerJob']
nSeed = config_dict['nSeed']
dtCalc = config_dict['dtCalc']
dtSave = config_dict['dtSave']
bufferLengthExport = config_dict['bufferLengthExport']
seedSpacing = config_dict['seedSpacing']
logLevel = config_dict['logLevel']
testSwitchHorizontal = config_dict['testSwitchHorizontal']
testSwitchVertical = config_dict['testSwitchVertical']

exportVariables = config_dict['exportVariables']
newVariables = config_dict['newVariables']
modelConfigDict = config_dict['modelConfigDict']

tracking_dir = tracking_dir_pre_str + outputDirName


# -----------------------------------------------
# With new input file format (yearly files), need new way of determing days in a year.  Use output of "ncdump" and regex capturing.
# -----------------------------------------------
# Bash integration because I forgot that this was a python script.  But maybe some of this is easier in bash... like, I need to use ncdump.
# -----------------------------------------------
test_regex='\\s*ocean_time = UNLIMITED ; \\/\\/ \\(([0-9]+) currently\\)'
#test_regex='ocean_time = UNLIMITED ; \/\/ \(([0-9]+) currently\)'
#test_regex='\tocean_time = UNLIMITED ; \/\/ \(([0-9]+) currently\)'
#test_regex='\s*ocean_time = UNLIMITED ; \/\/ \(([0-9]+) currently\)'
# -----------------------------------------------
# -----------------------------------------------


roms_run_name = inputDir.split('/')[-1]

roms_run_dir_list= os.listdir(path=inputDir)
roms_run_dir_list.sort()    

roms_run_dir_list = [os.path.join(inputDir,item) for item in roms_run_dir_list]
#roms_run_dir_list = [inputDir + "/" + item for item in roms_run_dir_list]

last_year = base_year + len(roms_run_dir_list)-1
#last_year = int(base_year) + len(roms_run_dir_list)-1


#details_str = ''
#if run_details[0] != None:
#    for item in run_details:
#        details_str = details_str + f'_{item}'
#if len(details_str) == 0:
#    experiment_dir = f'{roms_run_name}_{base_year}-{last_year}'
#else:
#    experiment_dir = f'{roms_run_name}_{details_str}_{base_year}-{last_year}'
#

experiment_dir = f'{roms_run_name}_{base_year}-{last_year}'

outputDir = os.path.join(baseOutputDir, tracking_dir, experiment_dir)

path = Path(outputDir)
path.mkdir(parents=True, exist_ok=True)

path = Path(outputDir + "/z_logs")
path.mkdir(parents=True, exist_ok=True)

path = Path(outputDir + "/z_config_files")
path.mkdir(parents=True, exist_ok=True)
config_dir = str(path)

day_nudge_run=nSeed*seedSpacing
day_nudge_job=day_nudge_run*numRunsPerJob


# With the new (yearly) file format, I'll just look at the first file in a year's directory (since I have 2 or 3 years in a year directory, to account for future years in a given tracking seeding),
# using ncdump -h and a few other piped commands to grab the line containing the "ocean_time" dimension, which is the number of days in a given year for Jerome's ROMS output

days_per_year_list=[]
for ii in range(len(roms_run_dir_list)):

    # Bash integration because I can't think
    cmd = f'''testFile='{roms_run_dir_list[ii]}'
    testLine="$(ncdump -h $testFile)"
    echo "$testLine"'''

    # Apparently this is "horrible" according to the SO author I copied
    bash_output = subprocess.run(cmd,shell=True,text=True,capture_output=True,check=True)
    bash_output = bash_output.stdout
    
    re_data = re.search(test_regex, bash_output)

    #print(ii)

    days_in_current_year = int(re_data.group(1))
   
    #pdb.set_trace()

    days_per_year_list.append(days_in_current_year)


cumulative_days_per_year_list=[]
total_days=0
for ii in range(len(days_per_year_list)):
    total_days = total_days + days_per_year_list[ii]
    cumulative_days_per_year_list.append(total_days)

#print(total_days)

#last_seed_day = total_days-int(modelConfigDict['drift:max_lifespan_days'])
last_seed_day = total_days-modelConfigDict['drift:max_lifespan_days']

#num_jobs=round((last_seed_day + (day_nudge_job-1))/day_nudge_job)
num_jobs = int(np.ceil(last_seed_day/day_nudge_job))  # why was I doing the previous thing?

dayNudge=0
runYear=0

#print(f"Particle Lifetime: {modelConfigDict['drift:max_lifespan_days']}")

for ii in range(num_jobs):

    dayNudge = ii*day_nudge_job

    if (dayNudge > cumulative_days_per_year_list[runYear]):
        for gg in range(runYear, len(roms_run_dir_list)):
            if (dayNudge > cumulative_days_per_year_list[gg]):
                runYear=((runYear+1))
            else:
                break

    startNudgeList=[]
    for jj in range(numRunsPerJob):
        currentNudge = dayNudge + day_nudge_run*jj
        if (currentNudge + day_nudge_run - seedSpacing) > last_seed_day:
            #print('FINAL NUDGE!')
            #print(f'last_seed_day = {last_seed_day}')
            break
        startNudgeList.append(currentNudge)

    if len(startNudgeList) == 0:
        print('startnudge breaking')
        break

    #print(startNudgeList)
    #print('lastNudge:')
    #print(startNudgeList[-1])
    #print('yearDays:')
    #print(cumulative_days_per_year_list[-2])
    #print(cumulative_days_per_year_list[-1])

    runYear0=0
    for cumulativeDays in cumulative_days_per_year_list:

        if startNudgeList[0] > cumulativeDays:
            runYear0=runYear0 + 1
        else:
            break

    runYear=runYear0
    
    #jobDir = roms_run_dir_list[runYear]

    #singleDirSwitchList=[]
    #jobDirList=[]
    #for nudge in startNudgeList:
    #    if runYear < len(roms_run_dir_list):
    #        if nudge > cumulative_days_per_year_list[runYear]:
    #            runYear = runYear+1
    #    if runYear+1 == len(roms_run_dir_list):
    #        singleDirSwitchList.append(1)
    #    else:
    #        singleDirSwitchList.append(0)
    #    jobDirList.append(roms_run_dir_list[runYear])


    logString="$(printf %02d ${dtCalc})_$(printf %04d ${dtSave})_$(printf %03d ${bufferLengthExport})_$(printf %02d ${nSeed})_$(printf %02d ${nRuns})_$(printf %06d ${ii})"
    
    config_file_name = "configFile_{c:03d}.config.yaml".format(c=ii)
    #config_file_name = "nRunsPerNode_{a:02d}_nSeed_{b:02d}_{c:03d}.config.yaml".format(a=numRunsPerJob, b=nSeed, c=ii)

    config_file = os.path.join(config_dir,config_file_name)

    shutil.copyfile(config_template_file, config_file)
            
    cd = {}

    cd["seedSpacing"] = seedSpacing
    cd["runCalc"] = dtCalc
    cd["runSave"] = dtSave
    cd["bufferLengthExport"] = bufferLengthExport
    cd["zznumberOfSeeds"]= nSeed 
    cd['numRunsPerJob'] = numRunsPerJob

    #cd["jobDir"] = jobDir
    #cd["zjobDirList"] = jobDirList
    #cd["dirListTotal"] = roms_run_dir_list 
    
    cd["zstartNudgeList"] = startNudgeList
    cd["inputDir"] = inputDir
    cd["outputDir"] = outputDir

    cd["baseYear"] = base_year

    cd['exportVariables'] = exportVariables
    cd['newVariables'] = newVariables
    cd['modelConfigDict'] = modelConfigDict

    cd['logLevel'] = logLevel
    cd['testSwitchVertical'] = testSwitchVertical
    cd['testSwitchHorizontal'] = testSwitchHorizontal



    with open(config_file, 'w') as outfile:
        yaml.dump(cd, outfile, default_flow_style=False)


    # If any second argument is provided to this script, only produce a single config file (for test runs)
    if args.testSwitch != None:
        break


