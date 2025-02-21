#!/bin/bash

inFileStr="wc15"

runDirStr="/home/blaughli/tracking_project_v2/misc/dummy_files/zero_vel_10days/single_run/Run_1995"
##runDirStr="/data/blaughli/jerome_projections_1990_2019/WC15N_GFDLTV/Run_"

catDir="/home/blaughli/concatenate_test"

catFile="wc15n_zero_vel_1995.nc"

ncrcat $runDirStr/$inFileStr* $catDir/$catFile
