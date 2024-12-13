#!/bin/bash

# Set the number of nodes to use for the job.  I usually check 2 things: call "sinfo" to see how many nodes are idle, and call "ls -l | wc -l" (subtract 1 from the result) in
# the z_config directory within the run directory, to see how many config files (ie jobs) there are for the run.  I'm typically generating 18 config files for 30 years of 
# ROMS output, so using 9 nodes means each node just runs 2 jobs, so the only way to complete the overall job more quickly would be to use 18 nodes, which I would probably
# only do if cluster activity seemed really low and/or after talking to lab mates

numNodes=9
#numNodes=6
#numNodes=2


runBaseDir=$1

callingDir="$(pwd)"

# This directory is populated by running the python script that generates the config files
runDirs=($runBaseDir/*)
    
extraArgs=""

numFiles=0

# Need to determine "serialSize" and "numNodesAtMaxSerial" before starting the main loop
for ii in "${!runDirs[@]}"
do
    runDir=${runDirs[$ii]}
    configFiles=($runDir/z_config_files/*)
    (( numFiles+=${#configFiles[@]} ))
done 


serialSize=$(( ($numFiles+$numNodes-1)/$numNodes ))
# Store the number of nodes we want run with $serialSize files (the rest we'll run at $serialSize-1 files)
numNodesAtMaxSerial=$(( $numNodes*(1-$serialSize)+$numFiles ))


counterRun=0
counterNode=0

for ii in "${!runDirs[@]}"
do

    runDir=${runDirs[$ii]}

    configFiles=($runDir/z_config_files/*)

    #numFiles=${#configFiles[@]}

    # rounding up
    #serialSize=$(( ($numFiles+$numNodes-1)/$numNodes ))

    # Store the number of files we want processed at the $serialSize (the rest we'll run at $serialSize-1)
    #numNodesAtMaxSerial=$(( $numNodes*(1-$serialSize)+$numFiles ))

#    counterRun=0
#    counterNode=0

#    extraArgs=""

    for jj in "${!configFiles[@]}"; 
    do

        (( counterRun ++ ))

        configFile=${configFiles[$jj]}
       
        configFileNum=$jj

        jobNum=$(sbatch --parsable --export="ALL,configFile=$configFile,callingDir=$callingDir,configFileNum=$configFileNum,runDir=$runDir" $extraArgs opendrift_caller.bash) 

        # For testing, maybe use extraArgs="--afterok", which will kill the whole job if something fails.  but that will help with time wasted monitoring.          
        # For production, use "--afterany", so that if a single job fails, can re-run later using the config file (that's part of the beauty of the config file approach)
        
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
done

