#!/bin/bash

# Making symbolic links for sample output, to test connectivity from a run that failed

#fileSourceDir="/home/blaughli/tracking_project_v2/x_old/tracking_to_do_full_nR_20_nS_15_kick_0p1/single_model_physicsOnly_1995-2020_kick_0p1"

fileSourceDir=$1

declare -a fileNameMasks

fileNameMasks[0]="tracking_output_configFile_000"
fileNameMasks[1]="tracking_output_configFile_001"
fileNameMasks[2]="tracking_output_configFile_002"
fileNameMasks[3]="tracking_output_configFile_003"
fileNameMasks[4]="tracking_output_configFile_004"
fileNameMasks[5]="tracking_output_configFile_005"
fileNameMasks[6]="tracking_output_configFile_006"
fileNameMasks[7]="tracking_output_configFile_007"
fileNameMasks[8]="tracking_output_configFile_008"
fileNameMasks[9]="tracking_output_configFile_009"

#symLinkDir="/home/blaughli/tracking_project_v2/tracking_to_do_sample/sample"

symLinkDir="/data03/blaughli/tracking_output/sample_10/sample"

mkdir -p $symLinkDir

for fileNameMask in ${fileNameMasks[@]}; do
    for entry in "$fileSourceDir"/"$fileNameMask"*; do
        fileName=$(echo "$entry" | awk -F/ '{print $NF}')
        ln -s "$entry" $symLinkDir/$fileName 
    done
done

