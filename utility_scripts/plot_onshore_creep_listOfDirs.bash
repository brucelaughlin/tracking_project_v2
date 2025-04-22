#!/bin/bash

dirNamePre=$1

fileNamePre=tracking_output_configFile_

for dir in ${PWD}/$dirNamePre*; do
    #echo "$dir"
    plot_onshore_creep $dir/$fileNamePre*
done
