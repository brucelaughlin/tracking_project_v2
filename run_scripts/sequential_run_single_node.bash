#!/bin/bash

# Run jobs sequentially on single machine

runDir=$1
configFiles=($runDir/z_config_files/*)
for configFile in "${configFiles[@]}"; 
do
    ./opendrift_caller.bash $configFile $runDir
done

