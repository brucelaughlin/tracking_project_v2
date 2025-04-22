#!/bin/bash

# Now trying to parallelize connectivityHist calculations.  My connectivityHist calculation is actually just a 2D histogram builder,
# with normalization happening later in a plotting script.  So, can just compute a bunch of histograms for various sets of 
# tracking files and then combine them.

# Note: Many thanks to The Paul for his counsel, which led to the switch towards a config-file appraoch.  The serialization
# described and implemented below was also his idea, and serialSize as a ceil calculation was his doing.

#runBaseDir="/data/blaughli/tracking_output/baseYear_1995_hindcast"

# Set niceLevel to 2000000000 to be as nice as possible.  0 to be greedy
niceLevel=0
#niceLevel=2000000000

runDir="$(realpath $1)"
#runDir=$1

connectivityHistJobFileDir="$runDir/y_connectivityHist_job_files"

mkdir -p "$connectivityHistJobFileDir"

# Set the number of nodes to use for the job (8 seems ok.. everyone else uses 8!)
#numNodes=9 # just works out that we have 18 jobs for now, so 9 splits them evenly between current and queued jobs
numNodes=$2

# Set the number of connectivityHist calculations to run on a node.  We have 48 cores per node on Atlantic, but I worry about using all of them...
# Also might get memory overflow, we'll see!
#numCoresUse=12
numCoresUse=44

# Location of THIS script (thanks SO)
callingDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# Need to determine "serialSize" and "numNodesAtMaxSerial" before starting the main loop
trackingFileList=($runDir/*.nc)
numTrackingFiles=0
(( numTrackingFiles+=${#trackingFileList[@]} ))

# Want to make "intermediate" job files which specify which tracking files a given job will process, so that I can use
# the bones of my original opendrift slurm submit scripts.

# Ceil division
(( numJobFiles = (numTrackingFiles+numCoresUse-1)/numCoresUse ))

jobFileNum=0

for trackingFileNum in $(seq 0 $numCoresUse $numTrackingFiles)
do
    jobFileName="$connectivityHistJobFileDir/jobFile_$(printf %02d $jobFileNum).txt"
    jobTrackingFileList=("${trackingFileList[@]:trackingFileNum:numCoresUse}")
    for trackingFile in "${jobTrackingFileList[@]}"
    do
        echo $trackingFile >> $jobFileName
    done
    ((jobFileNum++))
done


