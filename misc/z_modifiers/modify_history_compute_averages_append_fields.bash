#!/bin/bash

unmodifiedROMSDir=$(realpath "$1")
unmodifiedROMSArray=( "$unmodifiedROMSDir"/*nc )
modifiedROMSDir="${unmodifiedROMSDir}_AVERAGES"
mkdir -p "$modifiedROMSDir"
outputFile="${modifiedROMSDir}/u_v_avg.nc"
ncra -v u,v "${unmodifiedROMSDir}"/*.nc "$outputFile"
#ncra -v ubar,vbar "${unmodifiedROMSDir}"/*.nc "$outputFile"
sampleFile="${unmodifiedROMSArray[0]}"
ncks -A -v lon_rho,lat_rho,mask_rho "$sampleFile" "$outputFile"
