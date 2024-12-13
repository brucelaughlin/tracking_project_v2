#!/bin/bash

# Call the catch script on all tracking runs in a directory


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# Must specify directory containing run output directories to process
#--------------------------------------------------------------------------
# MUST CHANGE "baseYear" DEPENDING ON THE BASE YEAR OF THE TRACKING RUN
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#runBaseDir="/data/blaughli/copiedFiles/y_projections_1990_2020_KEEP/y_complete"
#baseYear=1990
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#runBaseDir="/data/blaughli/tracking_output/"
#baseYear=1988
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------


#baseYearDir=$1
#baseYearDirSplit=(${baseYearDir/_/ })
#baseYear=${baseYearDirSplit[1]}

baseYear=$1

runBaseDir="/data/blaughli/tracking_output/sanityTests/baseYear_$baseYear"
#runBaseDir="/data/blaughli/tracking_output/baseYear_$baseYear"

echo "$runBaseDir"
echo "$baseYear"

#callingDir="$(pwd)"
runDirArray=($runBaseDir/*)
   
for ii in "${!runDirArray[@]}"; do
    echo "${runDirArray[$ii]}"
    #screen -dmS "run_$ii" python catch_settled_points_rememberAllReleases_v2.py --trackingdir ${runDirArray[$ii]} --baseyear $baseYear
    python catch_settled_points_rememberAllReleases_v2.py --trackingdir ${runDirArray[$ii]} --baseyear $baseYear
done
wait

#wait # I don't know why I was using "wait" here.  I think the "&" above is required to make the python calls run in parallel

