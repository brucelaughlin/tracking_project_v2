
#------------------------------------------------------------------
pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_pld_20_30_z_two_file_test.p"
#pdf_file_name_pre = "pdf_data_output_releaseLoc_vs_settleTime_test3.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_test4_physics_only_AKs_1en5.p"
#------------------------------------------------------------------


#------------------------------------------------------------------
pdf_file_name = pdf_file_name_pre[0:-2] + "_swapped.p"
#------------------------------------------------------------------


import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import datetime

base_path = '/home/blaughli/tracking_project/'
pdf_raw_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/'
pdf_raw_file = pdf_raw_directory + pdf_file_name


file = open(pdf_raw_file,'rb')
bio_assigned_window,pdf_list_exposure_T,pdf_list_of_lists_O2,pdf_list_connectivity,pdf_list_settleTime,counter_array,box_num_mod,tick_positions,tick_labels,first_continent_box_num,oxygen_limit_list = pickle.load(file)
file.close()


pdf = pdf_list_settleTime[0]

pdf_colSum = np.sum(pdf,axis=0)


# something is weird in my assigned win

#bio_assigned_window = bio_assigned_window[1:]

plt.plot(pdf_colSum)
plt.plot(bio_assigned_window,c='r')
plt.show()

