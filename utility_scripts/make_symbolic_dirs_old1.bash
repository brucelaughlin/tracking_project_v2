#!/bin/bash

# Make the "symbolic link" directories, each with three links to three years except for the last two years.

fileSourceDir="/mnt/slough.pbsci.ucsc.edu/data/cedwards/WC15_hindcast_from_atlantic_202410_compressed"
#fileSourceDir="/data/cedwards/WC15_hindcast_from_atlantic_202410_compressed"
fileNameMask="wc15_fiechter_era5_glorys_avg_"

symLinkDir="/home/blaughli/symbolic_links_ROMS/WC15_hindcast_from_atlantic_202410_compressed/single_model"
#symLinkDir="/home/blaughli/symbolic_links_ROMS/WC15_hindcast_from_atlantic_202410_compressed"

yearFirst=1995
yearLast=2020
yearLastMinusOne=2019
yearLastMinusTwo=2018

for ii in $(seq $yearFirst $yearLast); do

    runYearPath=$symLinkDir/Run_$ii
        
    mkdir -p $runYearPath
    
    if (( ii == yearLast)); then 
        ln -s $fileSourceDir/$fileNameMask$ii.nc $runYearPath/$fileNameMask$ii.nc
    
    elif (( ii == yearLastMinusOne)); then 
        ln -s $fileSourceDir/$fileNameMask$ii.nc $runYearPath/$fileNameMask$ii.nc
        ln -s $fileSourceDir/$fileNameMask$((ii+1)).nc $runYearPath/$fileNameMask$((ii+1)).nc
    
    elif (( ii == yearLastMinusTwo)); then 
        ln -s $fileSourceDir/$fileNameMask$ii.nc $runYearPath/$fileNameMask$ii.nc
        ln -s $fileSourceDir/$fileNameMask$((ii+1)).nc $runYearPath/$fileNameMask$((ii+1)).nc
        ln -s $fileSourceDir/$fileNameMask$((ii+2)).nc $runYearPath/$fileNameMask$((ii+2)).nc

    else
        ln -s $fileSourceDir/$fileNameMask$ii.nc $runYearPath/$fileNameMask$ii.nc
        ln -s $fileSourceDir/$fileNameMask$((ii+1)).nc $runYearPath/$fileNameMask$((ii+1)).nc
        ln -s $fileSourceDir/$fileNameMask$((ii+2)).nc $runYearPath/$fileNameMask$((ii+2)).nc
        ln -s $fileSourceDir/$fileNameMask$((ii+3)).nc $runYearPath/$fileNameMask$((ii+3)).nc

    fi

done


