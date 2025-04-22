#!/bin/bash

# Make the "symbolic link" directories, each with three links to three years except for the last two years.
# Now with list of input directories.

#fileSourceDirList=("/home/cedwards/fiechter/WC15_era5_glorys_hindcast")
#fileSourceDirList=("${fileSourceDirList[@]}" "/data03/cae/fiechter/WC15_era5_glorys_hindcast")
fileSourceDirList=(/mnt/atlantic.pmc.ucsc.edu/home/cedwards/fiechter/WC15_era5_glorys_hindcast)
fileSourceDirList=("${fileSourceDirList[@]}" "/mnt/atlantic.pmc.ucsc.edu/data03/cae/fiechter/WC15_era5_glorys_hindcast")


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


