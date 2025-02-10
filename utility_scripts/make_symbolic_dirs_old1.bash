#!/bin/bash

# Make the "symbolic link" directories, each with three links to three years except for the last two years.

fileSourceDir_1="/home/cedwards/fiechter/WC15_era5_glorys_hindcast"
fileSourceDir_2="/data03/cae/fiechter/WC15_era5_glorys_hindcast"


#fileSourceDir="/mnt/slough.pbsci.ucsc.edu/data/cedwards/WC15_hindcast_from_atlantic_202410_compressed"

fileNameMask="wc15_fiechter_era5_glorys_avg_"

#symLinkDir="/home/blaughli/symbolic_links_ROMS/WC15_hindcast_from_atlantic_202410_compressed/single_model"
#symLinkDir="/home/blaughli/symbolic_links_ROMS/WC15_hindcast_from_atlantic_202410_compressed"
symLinkDir="/home/blaughli/symbolic_links_ROMS/WC15_hindcast_from_atlantic_202410_unpacked"

yearFirst=1995
yearLast=2020

mkdir - p $symLinkDir

for ii in $(seq $yearFirst $yearLast); do

    ln -s $fileSourceDir/$fileNameMask$ii.nc $symLinkDir/$fileNameMask$ii.nc

done


