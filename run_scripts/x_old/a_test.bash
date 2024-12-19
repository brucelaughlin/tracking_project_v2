#!/bin/bash

testFile='/data/cedwards/WC15_hindcast_from_atlantic_202410_compressed/wc15_fiechter_era5_glorys_avg_1995.nc'

testRegex='\s*ocean_time = UNLIMITED ; \/\/ \(([0-9]+) currently\)'
#testRegex='([0-9]+)'

testLine="$(ncdump -h $testFile | head -n 3 | tail -1)"

if [[ $testLine =~ $testRegex  ]] ; then
    num="${BASH_REMATCH[1]}"
fi

echo "$num"

