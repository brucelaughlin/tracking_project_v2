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

connectivityHistJobFileDir="$runDir/y_connectivity_hist_job_files"

mkdir -p "$connectivityHistJobFileDir"

# Set the number of nodes to use for the job (8 seems ok.. everyone else uses 8!)
#numNodes=9 # just works out that we have 18 jobs for now, so 9 splits them evenly between current and queued jobs
numNodes=$2

# Set the number of connectivityHist calculations to run on a node.  We have 24 cores per node on Atlantic...
# To check memory usage, run htop after ssh-ing into your allocated node
#numCoresUse=12
#numCoresUse=20 # Testing
#numCoresUse=24 # Testing
numCoresUse=22

# Location of THIS script (thanks SO)
callingDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# Need to determine "serialSize" and "numNodesAtMaxSerial" before starting the main loop
trackingFileList=($runDir/*.nc)
numTrackingFiles=0
(( numTrackingFiles+=${#trackingFileList[@]} ))

# Make "intermediate" job files which specify which tracking files a given job will process, so that I can use
# the bones of my original opendrift slurm submit scripts.

# Ceil division
(( numJobFiles = (numTrackingFiles+numCoresUse-1)/numCoresUse ))

jobFileNum=0

for trackingFileNum in $(seq 0 $numCoresUse $numTrackingFiles)
do
    jobFileName="$connectivityHistJobFileDir/jobFile_$(printf %03d $jobFileNum).txt"
    jobTrackingFileList=("${trackingFileList[@]:trackingFileNum:numCoresUse}")
    for trackingFile in "${jobTrackingFileList[@]}"
    do
        echo $trackingFile >> $jobFileName
    done
    ((jobFileNum++))
done

# Need to determine "serialSize" and "numNodesAtMaxSerial" before starting the main loop
connJobFileList=($connectivityHistJobFileDir/*)
numConnJobFiles=0
(( numConnJobFiles+=${#trackingFileList[@]} ))




# Say we have 8 config files (ie 8 slurm calls, ie 8 nodes' worth of work), and 3 nodes to work with.  We run in a serial way, so that
# file 2 starts as soon as file 1 finishes, file 5 starts as soon as file 4 finishes, etc.  Looks like this:
# 1 -> 2 -> 3
# 4 -> 5 -> 6
# 7 -> 8
# Note that Paul re-configured his version of this so that 1,2,3 are running together, and 1->4->7 is a serial chain (so that if
# a job is killed, completed files are in the expected order rather than being "shuffled" by the more primitive serial algorithm.

# <serialSize> is the maximum chain length in the above diagram.  It's equal to np.ceil(<numConnJobFiles>,<numNodes>), but we need to code
# our own version of "ceil" using bash's integer division, which is floor division (ie 9/10 = 0 in bash).
serialSize=$(( ($numConnJobFiles+$numNodes-1)/$numNodes ))
# Store the number of nodes that we'll use <serialSize> times (the rest we'll use <serialSize-1> times).  In the above diagram, 2/3 nodes
# get used <serialSize> (ie three) times, so <numNodesAtMaxSerial>=2 in this case.
# (I worked this math out on scratch paper and need to explain it better in future documentation)
numNodesAtMaxSerial=$(( $numNodes*(1-$serialSize)+$numConnJobFiles ))


extraArgs=""
counterRun=0
counterNode=0


for jj in "${!connJobFileList[@]}"; 
do
    (( counterRun ++ ))
    connJobFile=${connJobFileList[$jj]}
    connJobFileNum=$jj

    #echo "$connJobFile"
    #######/home/blaughli/tracking_project_v2/processing/settlement/y_parallel/sbatch_connHist_calc_call.bash $connJobFile


    jobNum=$(sbatch --nice=$niceLevel --parsable --export="ALL,connJobFile=$connJobFile,callingDir=$callingDir" $extraArgs /home/blaughli/tracking_project_v2/processing/settlement/y_parallel/sbatch_connHist_calc_call.bash)
    
#    # Test - run once
#    if (( $counterRun > 0 )); then                                                                                                                           
#        break
#    fi

    # For testing, maybe use extraArgs="--afterok", which will kill the whole job if something fails. 
    # For production, use "--afterany", so that if a single job fails, we can run it again using the associated config file 

    #extraArgs="-d afterany:$jobNum"                                                                                                                               
    extraArgs="-d afterok:$jobNum"                                                                                                                               
                                                                                                                                                                     
    if [[ $counterRun == $serialSize ]]; then                                                                                                                           
        counterRun=0                                                                                                                                                    
        extraArgs=""
        (( counterNode ++ )) 
        if [[ $counterNode == $numNodesAtMaxSerial ]]; then
            (( serialSize -- ))
        fi           
    fi                                                                                                                                                               
    
done

