#!/bin/bash

# Call the "modify" script on all settle algorithm aoutput in a directory


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# Must specify directory containing output directories to process
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
baseDir="/home/blaughli/tracking_project/practice/bounding_boxes/final_locations/z_output/z_pre_swap/z_swapped"
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------



dirArray=($baseDir/*.npz)
   
for f_file in "${dirArray[@]}"; do

    python plot_con_seasonal_pDrakeCompare_v1.py "$f_file"
    #python plot_con_seasonal_v9_patrick.py "$f_file"
    #python plot_con_seasonal_v8.py "$f_file"

    #break

done

