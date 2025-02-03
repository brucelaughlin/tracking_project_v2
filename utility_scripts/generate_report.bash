#!/bin/bash

runDir=$1
numNodes=$2

scriptPath="/home/blaughli/tracking_project_v2/misc/check_completion.py "

reportDir="x_report"
fileName="report.txt"
filePath=$runDir/$reportDir/$fileName
mkdir -p $runDir/$reportDir

python $scriptPath $runDir $numNodes > $filePath







