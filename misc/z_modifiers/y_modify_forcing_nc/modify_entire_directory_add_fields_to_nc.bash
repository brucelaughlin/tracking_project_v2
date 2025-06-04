#!/bin/bash

fileDir="$(realpath $1)"

fileList=($fileDir/*.nc)

for ii in "${!fileList[@]}";
do
    #echo "${fileList[$ii]}"
    ncks -A -v sea_floor_depth_below_sea_level,land_binary_mask /home/blaughli/symbolic_links_ROMS/z_testing/Mercator/global-reanalysis-phy-001-030-daily_1993_addVars.nc "${fileList[$ii]}"
done
