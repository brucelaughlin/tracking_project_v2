#!/bin/bash

#SBATCH --partition=memory-95GB
#SBATCH --job-name opendrift
#SBATCH --mail-type=ALL
#SBATCH --mail-user=blaughli@ucsc.edu

loggerLevel="DEBUG"

#configFile=$1
#runDir=$2

cd "$callingDir"
    
configFileTrim1="${configFile##*/}"
configFilePrefix="configFile_"
configFileTrim2=${configFileTrim1#"$configFilePrefix"}
configFileNum="${configFileTrim2/.config.yaml/}"

###slurmOutFile="slurm-configFile_$configFileNum.out"

#----------------------------------------------------------------------------------------------------------------
# So I think here I need to loop over the number of "job directories" specified in the config file
#----------------------------------------------------------------------------------------------------------------
# THIS FEELS SO HACK-Y
# So, our config file needs to always have the same order of dictionary keys, since my hacked solution below
# relies on it

# Hacked way of getting the number of opendrift calls we'll make on the current node.  Note that this approach
# REQUIRES that modifications to the config file preserve the order of 'zstartNudgeList' and 'zznumberOfSeeds'
# in the text of the yaml file
jobNudgeList=($(sed 's/[.-]//g' <<< $(awk '/zstartNudgeList/,/zznumberOfSeeds/' $configFile)))
jobNudgeList=("${jobNudgeList[@]:1:${#jobNudgeList[@]}-3}")
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

for jobRunNum in "${!jobNudgeList[@]}"; do

    jobRunString="${configFilePrefix}${configFileNum}_job_$(printf %02d ${jobRunNum})"
    logFilePre="$jobRunString.driftlog"
    logFile="${runDir}/z_logs/$(basename $logFilePre)"

    echo "$(hostname)" > "$logFile"

    python /home/blaughli/tracking_project_v2/run_scripts/opendrift_run.py --configfile $configFile --jobrunnumber $jobRunNum --level $loggerLevel --jobrunstring $jobRunString &>> "$logFile" &
    

done
wait # I don't know why I was using "wait" here.  I think the "&" above is required to make the python calls run in parallel

