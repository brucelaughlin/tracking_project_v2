#!/bin/bash

#cd "$callingDir"

loggerLevel="DEBUG"

configFile=$1
#configFileNum=$2
runDir=$2
#runDir=$3


#----------------------------------------------------------------------------------------------------------------
# So I think here I need to loop over the number of "job directories" specified in the config file
#----------------------------------------------------------------------------------------------------------------
# THIS FEELS SO HACK-Y
# So, our config file needs to always have the same order of dictionary keys, since my hacked solution below
# relies on it

jobNudgeList=($(sed 's/[.-]//g' <<< $(awk '/zstartNudgeList/,/zznumberOfSeeds/' $configFile)))
jobNudgeList=("${jobNudgeList[@]:1:${#jobNudgeList[@]}-3}")
#jobDirList=($(sed 's/[.-]//g' <<< $(awk '/zjobDirList/,/zznumberOfSeeds/' $configFile)))
#jobDirList=("${jobDirList[@]:1:${#jobDirList[@]}-3}")
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

for jobRunNum in "${!jobNudgeList[@]}"; do

    configFileTrim1="${configFile##*/}"
    configFilePrefix="configFile_"
    configFileTrim2=${configFileTrim1#"$configFilePrefix"}
    configFileNum="${configFileTrim2/.config.yaml/}"

    #jobRunString=printf "%s%d_%s" "$configFilePrefix" "$configFileNum" "$(printf %02d ${jobRunNum})"
    jobRunString="${configFilePrefix}${configFileNum}_job_$(printf %02d ${jobRunNum})"
    #jobRunString="${configFile/.config.yaml/}_$(printf %02d ${jobRunNum})"
    logFilePre="$jobRunString.driftlog"
    #logFilePre="${configFile/.config.yaml/}_configFileNum_$(printf %02d ${configFileNum})_configFileJob_$(printf %02d ${jobRunNum}).driftlog"
    logFile="${runDir}/z_logs/$(basename $logFilePre)"

    #echo "$logFile"
    echo "$(hostname)" > "$logFile"

    python opendrift_run.py --configfile $configFile --jobrunnumber $jobRunNum --level $loggerLevel --jobrunstring $jobRunString &>> "$logFile" &
    
    ###python opendrift_run.py --configfile $configFile --jobrunnumber $jobRunNum --level $loggerLevel &>> "$logFile" &

done
wait # I don't know why I was using "wait" here.  I think the "&" above is required to make the python calls run in parallel

