#!/bin/bash

# Make the "symbolic link" directories, each with three links to three years except for the last two years.

#numSourceDir=2

#fileSourceDir1="/home/cedwards/fiechter/WC15_era5_glorys_hindcast"
#fileSourceDir2="/data03/cae/fiechter/WC15_era5_glorys_hindcast"

#fileSourceDirList=("$fileSourceDir1")
fileSourceDirList=("/home/cedwards/fiechter/WC15_era5_glorys_hindcast")
fileSourceDirList=("${fileSourceDirList[@]}" "/data03/cae/fiechter/WC15_era5_glorys_hindcast")


#declare -a fileSourceDirList=()
#for ii in $(seq 1 $numSourceDir) 
#do
#    fileSourceDirList=("${fileSourceDirList[@]}" "")

#fileSourceDir="/mnt/slough.pbsci.ucsc.edu/data/cedwards/WC15_hindcast_from_atlantic_202410_compressed"





#fileNameMask="wc15_fiechter_era5_glorys_avg_"

#symLinkDir="/home/blaughli/symbolic_links_ROMS/WC15_hindcast_from_atlantic_202410_compressed/single_model"
#symLinkDir="/home/blaughli/symbolic_links_ROMS/WC15_hindcast_from_atlantic_202410_compressed"
symLinkDir="/home/blaughli/symbolic_links_ROMS/WC15_hindcast_from_atlantic_202410_unpacked"

#yearFirst=1995
#yearLast=2020

mkdir -p $symLinkDir

#fileList1='ls $fileSourceDir1/*.nc'
#fileList2='ls $fileSourceDir2/*.nc'

for fileSourceDir in "${fileSourceDirList[@]}"
do
    #echo "$fileSourceDir"

    #fileList='ls $fileSourceDir/*.nc'
    #for fileName in "${fileList1[@]}"
    for fileName in "$fileSourceDir"/*.nc
    do 
        #echo "$fileName"
        ln -s $fileName $symLinkDir/"$(basename $fileName)"
    done
done






#for ii in $(seq $yearFirst $yearLast); do
#    ln -s $fileSourceDir/$fileNameMask$ii.nc $symLinkDir/$fileNameMask$ii.nc
#done


