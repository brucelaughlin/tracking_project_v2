#!/bin/bash

# Making symbolic links for sample output, to test connectivity from a run that failed

#fileSourceDir="/home/blaughli/tracking_project_v2/x_old/tracking_to_do_full_nR_20_nS_15_kick_0p1/single_model_physicsOnly_1995-2020_kick_0p1"

fileSourceDir=$1

fileNameMask="tracking_output_configFile_000"

#symLinkDir="/home/blaughli/tracking_project_v2/tracking_to_do_sample/sample"

symLinkDir="/data03/blaughli/tracking_output/sample/sample"

mkdir -p $symLinkDir

for entry in "$fileSourceDir"/"$fileNameMask"*; do
    
    fileName=$(echo "$entry" | awk -F/ '{print $NF}')

    ln -s "$entry" $symLinkDir/$fileName 

done


