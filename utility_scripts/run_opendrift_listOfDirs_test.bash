#!/bin/bash

dirNamePre=$1
numNodesPerJob=$2

for dir in ${PWD}/$dirNamePre*; do
#for dir in ${listOfDirs[@]}; do
    tests_submit_slurm_jobs_opendrift $dir $numNodesPerJob
    #echo "$dir"
done
