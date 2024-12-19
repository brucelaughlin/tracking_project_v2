#!/bin/bash

#cd "$callingDir"

#----------------------------------------------------------------------------------------------------------------
# So I think here I need to loop over the number of "job directories" specified in the config file
#----------------------------------------------------------------------------------------------------------------
# THIS FEELS SO HACK-Y
# So, our config file needs to always have the same order of dictionary keys, since my hacked solution below
# relies on it

configFile="/home/blaughli/tracking_project_v2/tracking_to_do_1990_readerBugTest/WC15N_GFDLTV_physicsOnly_1990-2019/z_config_files/nRunsPerNode_15_nSeed_20_001.config.yaml"

jobDirList=($(sed 's/[.-]//g' <<< $(awk '/zjobDirList/,/zznumberOfSeeds/' $configFile)))
jobDirList=("${jobDirList[@]:1:${#jobDirList[@]}-3}")
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

for jobRunNum in "${!jobDirList[@]}"; do

    #logFilePre="${configFile/.config.yaml/}_configFileNum_$(printf %02d ${configFileNum})_configFileJob_$(printf %02d ${jobRunNum}).driftlog"
    #logFile="${runDir}/z_logs/$(basename $logFilePre)"

    #echo "$(hostname)" > "$logFile"

    echo "${jobDirList[$jobRunNum]}"

done
wait # I don't know why I was using "wait" here.  I think the "&" above is required to make the python calls run in parallel

