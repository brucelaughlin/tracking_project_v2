#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_oneFileTest_swapped.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_tenFileTest_swapped.p"
#pdf_file_name_pre = "pdf_data_output_releaseLoc_vs_settleTime_test3_swapped.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_run_test3_swapped.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_test4_physics_only_AKs_1en5_swapped.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_ranges_O2_pH___pld_20_30_z_two_file_test_swapped.npz"
pdf_file_name_pre = "pdf_data_output_seasonal_ranges_O2_pH___pld_20_29_z_two_file_test_swapped.npz"
#pdf_file_name_pre = "pdf_data_output_seasonal_ranges_O2_pH___pld_90_149_drift_150_physics_only_AKs_1en5_v1_swapped.npz"
#------------------------------------------------------------------

# Testing - choose "test_dex", corresponding to entries from 
###pH_limit_list = [7.6,7.7,7.8,8,8.5,9]
#pH_limit_list = [7.5,7.6,7.7,7.8,7.9,8,8.1,8.2,8.3]

#test_dex = 2
test_dex = 5
#test_dex = 8
#------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import datetime

#------------------------------------------------------------------
#pdf_file_name = pdf_file_name_pre[0:-2] + "_swapped.p"
pdf_file_name = pdf_file_name_pre

base_path = '/home/blaughli/tracking_project/'
pdf_raw_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/'

pdf_raw_file = pdf_raw_directory + pdf_file_name
#------------------------------------------------------------------

d = np.load(pdf_raw_file)

pdf_arrays_T = d['pdf_arrays_T']
pdf_arrays_O2 = d['pdf_arrays_O2']
pdf_arrays_pH = d['pdf_arrays_pH']
pdf_arrays_connectivity = d['pdf_arrays_connectivity']
pdf_arrays_settleTime = d['pdf_arrays_settleTime']
#counter_array = d['counter_array']
#oxygen_limit_list = d['oxygen_limit_list']
pH_limit_list = d['pH_limit_list']
box_num_mod = d['box_num_mod']
tick_positions = d['tick_positions']
tick_labels = d['tick_labels']
first_continent_box_dex = d['first_continent_box_num']

#a = pdf_arrays_settleTime[0,:,:]
a = pdf_arrays_pH[test_dex,0,:,:]

plt.pcolormesh(a.T)
plt.show()

