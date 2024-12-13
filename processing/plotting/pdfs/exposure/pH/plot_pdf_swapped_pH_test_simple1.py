# v1: copied from v4 of oxygen


#------------------------------------------------------------------
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_oneFileTest_swapped.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_tenFileTest_swapped.p"
#pdf_file_name_pre = "pdf_data_output_releaseLoc_vs_settleTime_test3_swapped.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_run_test3_swapped.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_test4_physics_only_AKs_1en5_swapped.p"
pdf_file_name_pre = "pdf_data_output_seasonal_ranges_O2_pH___pld_20_30_z_two_file_test_swapped.npz"
#------------------------------------------------------------------

# -------------------------------------------------------------------
# Testing - choose 0-3 for test dex, corresponding to entries from 
###pH_limit_list = [7.6,7.7,7.8,8,8.5,9]
#pH_limit_list = [7.5,7.6,7.7,7.8,7.9,8,8.1,8.2,8.3]

test_dex = 8
#test_dex = 5
#------------------------------------------------------------------





#pdf_file_name = pdf_file_name_pre[0:-2] + "_swapped.p"
pdf_file_name = pdf_file_name_pre
#------------------------------------------------------------------

# Titles, labels, random params
#------------------------------------------------------------------
fig_param_title = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-90 day PLD\n"
#fig_param_title = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD\n"
#fig_param_title = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"
#fig_main_title_pre = "PDFs of cumulative number of days larvae were exposed to DO2 concentrations below {} mg/L (y-axis) vs Release Location (x-axis).\nGrouped according to season of release."






y_label_pre = "exposure time (days)"
#y_label_pre = "exposure time to [DO2] below {} mg/L (days)"

day_0_title_pre = "probability of 0 days of exposure"
#day_0_title_pre = "probability of 0 days of exposure to [DO2] below {} mg/L"

day_1_60_title_pre = "probability of 1-60 days of exposure  (settler release months:"
#day_1_60_title_pre = "probability of 1-60 days of exposure  (release months:"
#day_1_60_title_pre = "probability of 1-60 days of exposure"
#day_1_60_title_pre = "probabilities of 1-60 days of exposure"
#day_1_60_title_pre = "probabilities of 1-60 days of exposure to [DO2] below {} mg/L"

cbar_label = "probability"
#cbar_label = "Log base 10 value of probability"
cbar_fontSize = 15
cbar_nBins = 25


#------------------------------------------------------------------


import pickle
import numpy as np
import matplotlib.pyplot as plt
import datetime

#------------------------------------------------------------------
base_path = '/home/blaughli/tracking_project/'
pdf_raw_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/'

pdf_raw_file = pdf_raw_directory + pdf_file_name
#------------------------------------------------------------------

d = np.load(pdf_raw_file)

#pdf_list_exposure_T = d['pdf_list_exposure_T']
#pdf_list_of_lists_O2 = d['pdf_list_of_lists_O2']
pdf_list_of_lists_pH = d['pdf_list_of_lists_pH']
#pdf_list_connectivity = d['pdf_list_connectivity']
#pdf_list_settleTime = d['pdf_list_settleTime']
#counter_array = d['counter_array']
#oxygen_limit_list = d['oxygen_limit_list']
pH_limit_list = d['pH_limit_list']
box_num_mod = d['box_num_mod']
tick_positions = d['tick_positions']
tick_labels = d['tick_labels']
first_continent_box_dex = d['first_continent_box_num']

a = pdf_list_of_lists_pH[test_dex,0,:,:]

plt.pcolormesh(a.T)
plt.show()

