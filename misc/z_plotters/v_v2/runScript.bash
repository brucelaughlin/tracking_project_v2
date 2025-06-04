#!/bin/bash


#gridDir="/home/blaughli/tracking_project_v2/grid_data/mercator_diy_grid.npz"
gridDir="/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_noModification.npz"

floatDex=0

runDir="/data03/blaughli/tracking_output/Mercator_coastalCells_10dayTest_kick_0p0_defaultMask___onshore_uo/tracking_output_configFile_000_job_00.nc"
python plot_single_trajectory_Mercator_update.py "$runDir" "$gridDir" "$floatDex" &

runDir="/data03/blaughli/tracking_output/Mercator_coastalCells_10dayTest_kick_0p0_customMask___onshore_uo/tracking_output_configFile_000_job_00.nc"
python plot_single_trajectory_Mercator_update.py "$runDir" "$gridDir" "$floatDex" &
