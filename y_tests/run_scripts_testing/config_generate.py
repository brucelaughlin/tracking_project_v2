# Generate config files for Opendrift Runs

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

config_template_file = '/home/blaughli/tracking_project_v2/config_files_for_setup/z_templates/config_template.config.yaml'

parser = argparse.ArgumentParser()
parser.add_argument('configfile', type=str)
#parser.add_argument('--configfile', type=str)
args = parser.parse_args()

config_for_config_file = args.configfile

stream = open(config_for_config_file,'r')
#config_dict = yaml.load(stream, Loader=yaml.BaseLoader)
config_dict = yaml.safe_load(stream)
stream.close()

his_file_name_pre = config_dict['hisFileNamePre']
base_year = config_dict['baseYear']
#base_year = int(config_dict['baseYear'])
baseInputDirPre = config_dict['baseInputDirPre']
romsRunDirectories = config_dict['romsRunDirectories']
baseOutputDir = config_dict['baseOutputDir']
behavior_list = config_dict['behaviorList']
run_details = config_dict['runDetails']

numRunsPerJob = config_dict['numRunsPerJob']
nSeed = config_dict['nSeed']
dtCalc = config_dict['dtCalc']
dtSave = config_dict['dtSave']
bufferLength = config_dict['bufferLength']
seedSpacing = config_dict['seedSpacing']
logLevel = config_dict['logLevel']
singleParticleTestSwitch = config_dict['singleParticleTestSwitch']

#numRunsPerJob = int(config_dict['numRunsPerJob'])
#nSeed = int(config_dict['nSeed'])
#dtCalc = int(config_dict['dtCalc'])
#dtSave = int(config_dict['dtSave'])
#bufferLength = int(config_dict['bufferLength'])
#seedSpacing = int(config_dict['seedSpacing'])
#logLevel = int(config_dict['logLevel'])

exportVariables = config_dict['exportVariables']
newVariables = config_dict['newVariables']
modelConfigDict = config_dict['modelConfigDict']

for behavior in behavior_list:

    for roms_run_number in range(len(romsRunDirectories)):

        roms_run_name = romsRunDirectories[roms_run_number]

        base_input_dir = os.path.join(baseInputDirPre,roms_run_name)
        #base_input_dir = baseInputDirPre + roms_run_name

        roms_run_dir_list= os.listdir(path=base_input_dir)
        roms_run_dir_list.sort()    

        roms_run_dir_list = [os.path.join(base_input_dir,item) for item in roms_run_dir_list]
        #roms_run_dir_list = [base_input_dir + "/" + item for item in roms_run_dir_list]

        last_year = base_year + len(roms_run_dir_list)-1
        #last_year = int(base_year) + len(roms_run_dir_list)-1

        experiment_dir = f'{roms_run_name}_{behavior}_{base_year}-{last_year}'
        #experiment_dir = f'{roms_run_name}_{behavior}_{base_year}-{base_year+len(roms_run_dir_list)-1}'
        
        if run_details[0] != None:
        #if run_details[0] != '':
            for item in run_details:
                experiment_dir = experiment_dir + f'_{item}'

        outputDir = os.path.join(baseOutputDir,experiment_dir)

        path = Path(outputDir)
        path.mkdir(parents=True, exist_ok=True)

        path = Path(outputDir + "/z_logs")
        path.mkdir(parents=True, exist_ok=True)

        path = Path(outputDir + "/z_config_files")
        path.mkdir(parents=True, exist_ok=True)
        config_dir = str(path)

        day_nudge_run=nSeed*seedSpacing
        day_nudge_job=day_nudge_run*numRunsPerJob

        days_per_year_list=[]
        for ii in range(len(roms_run_dir_list)):
            run_year_list= os.listdir(path=roms_run_dir_list[ii])
            run_year_list.sort()
            run_year_list = [item for item in run_year_list if item[0:len(his_file_name_pre)] == his_file_name_pre]
            days_per_year_list.append(len(run_year_list))


        cumulative_days_per_year_list=[]
        total_days=0
        for ii in range(len(days_per_year_list)):
            total_days = total_days + days_per_year_list[ii]
            cumulative_days_per_year_list.append(total_days)

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
                print('breaking!')
                break

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

            singleDirSwitchList=[]
            jobDirList=[]
            for nudge in startNudgeList:
                if runYear < len(roms_run_dir_list):
                    if nudge > cumulative_days_per_year_list[runYear]:
                        runYear = runYear+1
                if runYear+1 == len(roms_run_dir_list):
                    singleDirSwitchList.append(1)
                else:
                    singleDirSwitchList.append(0)
                jobDirList.append(roms_run_dir_list[runYear])


            logString="$(printf %02d ${dtCalc})_$(printf %04d ${dtSave})_$(printf %03d ${bufferLength})_$(printf %02d ${nSeed})_$(printf %02d ${nRuns})_$(printf %06d ${ii})"
            
            config_file_name = "nRunsPerNode_{a:02d}_nSeed_{b:02d}_{c:03d}.config.yaml".format(a=numRunsPerJob, b=nSeed, c=ii)

            config_file = os.path.join(config_dir,config_file_name)

            shutil.copyfile(config_template_file, config_file)
                    
            cd = {}

            cd["seedSpacing"] = seedSpacing
            cd["runCalc"] = dtCalc
            cd["runSave"] = dtSave
            cd["bufferLength"] = bufferLength
            cd["zznumberOfSeeds"]= nSeed 

            cd["zjobDirList"] = jobDirList
            cd["dirListTotal"] = roms_run_dir_list 
            
            cd["startNudgeList"] = startNudgeList
            cd["outputDir"] = outputDir

            cd["behavior"] = behavior
            cd["baseYear"] = base_year

            cd['exportVariables'] = exportVariables
            cd['newVariables'] = newVariables
            cd['modelConfigDict'] = modelConfigDict

            cd['logLevel'] = logLevel
            cd['singleParticleTestSwitch'] = singleParticleTestSwitch



            with open(config_file, 'w') as outfile:
                yaml.dump(cd, outfile, default_flow_style=False)



