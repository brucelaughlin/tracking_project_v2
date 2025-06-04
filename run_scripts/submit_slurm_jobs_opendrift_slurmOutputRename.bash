#!/bin/bash

# Note: Many thanks to The Paul for his counsel, which led to the switch towards a config-file appraoch.  The serialization
# described and implemented below was also his idea, and serialSize as a ceil calculation was his doing.

#runBaseDir="/data/blaughli/tracking_output/baseYear_1995_hindcast"

niceLevel=2000000000

runDir="$(realpath $1)"
#runDir=$1

# Set the number of nodes to use for the job (8 seems ok.. everyone else uses 8!)
#numNodes=9 # just works out that we have 18 jobs for now, so 9 splits them evenly between current and queued jobs
numNodes=$2

# Location of THIS script (thanks SO)
callingDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# Need to determine "serialSize" and "numNodesAtMaxSerial" before starting the main loop
configFileList=($runDir/z_config_files/*)
numConfigFiles=0
(( numConfigFiles+=${#configFileList[@]} ))


# Say we have 8 config files (ie 8 slurm calls, ie 8 nodes' worth of work), and 3 nodes to work with.  We run in a serial way, so that
# file 2 starts as soon as file 1 finishes, file 5 starts as soon as file 4 finishes, etc.  Looks like this:
# 1 -> 2 -> 3
# 4 -> 5 -> 6
# 7 -> 8
# Note that Paul re-configured his version of this so that 1,2,3 are running together, and 1->4->7 is a serial chain (so that if
# a job is killed, completed files are in the expected order rather than being "shuffled" by the more primitive serial algorithm.

# <serialSize> is the maximum chain length in the above diagram.  It's equal to np.ceil(<numConfigFiles>,<numNodes>), but we need to code
# our own version of "ceil" using bash's integer division, which is floor division (ie 9/10 = 0 in bash).
serialSize=$(( ($numConfigFiles+$numNodes-1)/$numNodes ))
# Store the number of nodes that we'll use <serialSize> times (the rest we'll use <serialSize-1> times).  In the above diagram, 2/3 nodes
# get used <serialSize> (ie three) times, so <numNodesAtMaxSerial>=2 in this case.
# (I worked this math out on scratch paper and need to explain it better in future documentation)
numNodesAtMaxSerial=$(( $numNodes*(1-$serialSize)+$numConfigFiles ))


extraArgs=""
counterRun=0
counterNode=0


for jj in "${!configFileList[@]}"; 
do
    (( counterRun ++ ))
    configFile=${configFileList[$jj]}
    configFileNum=$jj

    jobNum=$(sbatch --nice=$niceLevel --parsable --export="ALL,configFile=$configFile,callingDir=$callingDir,runDir=$runDir" $extraArgs /home/blaughli/tracking_project_v2/run_scripts/sbatch_opendrift_call.bash) 
   
    configFileTrim1="${configFile##*/}"
    configFilePrefix="configFile_"
    configFileTrim2=${configFileTrim1#"$configFilePrefix"}
    configFileNum="${configFileTrim2/.config.yaml/}"

    # Copy and rename the slurm-##.out file so that we know which config file was used, in case of failure 
    ln -f "slurm-$jobNum.out" "slurm-$jobNum-configFileNum-$configFileNum.out" 
    
    # For testing, maybe use extraArgs="--afterok", which will kill the whole job if something fails. 
    # For production, use "--afterany", so that if a single job fails, we can run it again using the associated config file 

    extraArgs="-d afterany:$jobNum"                                                                                                                               
    #extraArgs="-d afterok:$jobNum"                                                                                                                               
                                                                                                                                                                     
    if [[ $counterRun == $serialSize ]]; then                                                                                                                           
        counterRun=0                                                                                                                                                    
        extraArgs=""
        (( counterNode ++ )) 
        if [[ $counterNode == $numNodesAtMaxSerial ]]; then
            (( serialSize -- ))
        fi           
    fi                                                                                                                                                               
    
done

