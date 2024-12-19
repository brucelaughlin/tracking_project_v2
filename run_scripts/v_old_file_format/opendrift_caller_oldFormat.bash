#!/bin/bash

#cd "$callingDir"

loggerLevel="DEBUG"

configFile=$1
configFileNum=$2
runDir=$3


#----------------------------------------------------------------------------------------------------------------
# So I think here I need to loop over the number of "job directories" specified in the config file
#----------------------------------------------------------------------------------------------------------------
# THIS FEELS SO HACK-Y
# So, our config file needs to always have the same order of dictionary keys, since my hacked solution below
# relies on it

jobDirList=($(sed 's/[.-]//g' <<< $(awk '/zjobDirList/,/zznumberOfSeeds/' $configFile)))
jobDirList=("${jobDirList[@]:1:${#jobDirList[@]}-3}")
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

for jobRunNum in "${!jobDirList[@]}"; do

    logFilePre="${configFile/.config.yaml/}_configFileNum_$(printf %02d ${configFileNum})_configFileJob_$(printf %02d ${jobRunNum}).driftlog"
    logFile="${runDir}/z_logs/$(basename $logFilePre)"

    echo "$(hostname)" > "$logFile"

    python opendrift_run_v2.py --configfile $configFile --jobrunnumber $jobRunNum --level $loggerLevel &>> "$logFile" &
    #python opendrift_run_v2.py --configfile $configFile --jobrunnumber $jobRunNum &>> "$logFile" &

done
wait # I don't know why I was using "wait" here.  I think the "&" above is required to make the python calls run in parallel

