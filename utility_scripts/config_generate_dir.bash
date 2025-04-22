#!/bin/bash

#configDir=$realpath $1

for file in $PWD/*; do
#for file in $configDir/*; do

    if test -f "$file"
    then
        python /home/blaughli/tracking_project_v2/run_scripts/config_generate.py $file
    fi

done
