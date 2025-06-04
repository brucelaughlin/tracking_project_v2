#!/bin/bash

#SBATCH --job-name connCalc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=blaughli@ucsc.edu

loggerLevel="DEBUG"

#connJobFile=$1

cd "$callingDir"
    
connJobFileTrim1="${connJobFile##*/}"
connJobFilePrefix="jobFile_"
connJobFileTrim2=${connJobFileTrim1#"$connJobFilePrefix"}
connJobFileNum="${connJobFileTrim2/.txt/}"


connDataDir="$(dirname "$connJobFile")"
logDir="$connDataDir/l_logs"


# --------
# temporary, for testing
logFilePath="/home/blaughli/aa_log.txt"
# --------

IFS=

while read -r line; do
    echo "File: $line" >> "$logFilePath"

    trackingFileTrim1="${line##*/}"
    trackingFileTrim2="${trackingFileTrim1/.nc/}"

    logFilePre="${trackingFileTrim2}.log"
    logFilePath="$logDir/$logFilePre"


    python /home/blaughli/tracking_project_v2/processing/settlement/y_parallel/connectivity_histogram_calc_specify_polygons.py --trackingfile $line --pdrakefileswitch 0 --polygonfile $polygonFile  &>> "$logFilePath" &
    #python /home/blaughli/tracking_project_v2/processing/settlement/y_parallel/a_conn_test.py --trackingfile $line --pdrakefileswitch 0 --pdrakepolygonswitch 1  &>> "$logFilePath" &
done < "$connJobFile"

    

wait # I don't know why I was using "wait" here.  I think the "&" above is required to make the python calls run in parallel

