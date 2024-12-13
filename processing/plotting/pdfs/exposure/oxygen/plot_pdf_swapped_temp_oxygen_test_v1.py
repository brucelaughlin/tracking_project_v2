# Copied from ""

#------------------------------------------------------------------
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_oneFileTest_swapped.p"
pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_tenFileTest_swapped.p"
#pdf_file_name_pre = "pdf_data_output_releaseLoc_vs_settleTime_test3_swapped.p"

#------------------------------------------------------------------

#------------------------------------------------------------------
#pdf_file_name = pdf_file_name_pre[0:-2] + "_swapped.p"
pdf_file_name = pdf_file_name_pre
#------------------------------------------------------------------


#------------------------------------------------------------------
fig_param_title = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"
fig_main_title_pre = "PDFs of cumulative number of days larvae were exposed to DO2 concentrations below {} mg/L (y-axis) vs Release Location (x-axis).\nGrouped according to season of release."
#fig_main_title_pre = "PDFs of cumulative number of days exposed to DO2 levels below {} mg/L (y-axis)\n vs Release Location (x-axis).\nGrouped according to season of release."

y_label_pre = "exposure time to DO2 below {} mg/L (days)"



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


file = open(pdf_raw_file,'rb')

pdf_list_exposure_T_source_swapped,pdf_list_of_lists_O2_source_swapped,pdf_list_connectivity_swapped,pdf_list_settleTime_swapped,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels,first_continent_box_dex,oxygen_limit_list = pickle.load(file)
#pdf_list_exposure_T_source_swapped,pdf_list_of_lists_O2_source_swapped,pdf_list_connectivity_swapped,pdf_list_settleTime_swapped,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels,first_continent_box_dex = pickle.load(file)
file.close()


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Testing - choose 0-3 for test dex, corresponding to entries from 
#oxygen_limit_list = [2.2,3.1,4.1,6] # REMOVE THIS AFTER RUNNING CATCH SCRIPT AGAIN
test_dex = 3
# -------------------------------------------------------------------
pdf_list_exposure_oxygen_source = pdf_list_of_lists_O2_source_swapped[test_dex]
# -------------------------------------------------------------------
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Dynamic title and labels
fig_main_title = fig_main_title_pre.format(oxygen_limit_list[test_dex])
y_label = y_label_pre.format(oxygen_limit_list[test_dex])
# -------------------------------------------------------------------
fig_fullTitle = fig_main_title + "\n" + fig_param_title
# -------------------------------------------------------------------




pdf_list_full = []
pdf_max_val = -999999
pdf_min_val = 999999

# -------------------------------------------------------------------
for pdf in pdf_list_exposure_oxygen_source[1:]:
    pdf_full = pdf.copy()


    # -------------------------------------------------------------------
    # Convert to probability
    row_sums = pdf_full.sum(axis=1)
    pdf_full = pdf_full / row_sums[:, np.newaxis]
    pdf_full = np.log10(pdf_full)
    # -------------------------------------------------------------------
    pdf_list_full.append(pdf_full)
    if np.amax(pdf_full) > pdf_max_val:
        pdf_max_val = np.amax(pdf_full)
    if np.amin(np.ma.masked_invalid(pdf_full)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf_full))
# -------------------------------------------------------------------


# Determined elsewhere (see/run "check_box_numbers.py")
#first_continent_box_dex = 17
#first_continent_box_dex = 20
num_dummy_lines = 1

n_rows = int(np.shape(pdf_list_exposure_oxygen_source[0])[0]) + num_dummy_lines
n_columns = int(np.shape(pdf_list_exposure_oxygen_source[0])[1])
X = np.arange(-0.5, n_rows, 1)
Y = np.arange(-0.5, n_columns, 1)

v_scale = 6



fig,axs = plt.subplots(2,2)

label_fontSize = 10
#label_fontSize = 8

for ii in range(len(pdf_list_full)):

    # 2D full plot (minus day 1)
    pdf_plot = pdf_list_full[ii]
    pdf_separated = np.empty((np.shape(pdf_plot)[0] + num_dummy_lines,np.shape(pdf_plot)[1]))
    pdf_separated[:] = np.nan
    pdf_separated[0:first_continent_box_dex,:] = pdf_plot[0:first_continent_box_dex,:]
    pdf_separated[first_continent_box_dex + num_dummy_lines:,:] = pdf_plot[first_continent_box_dex:,:]

    if ii == 1:
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        pname = "Winter (DJF)"
        #axs[0,0].title.set_text("Winter (DJF)")
        #axs[0,0].title.set_text("Winter (Dec,Jan,Feb)")
        axs[0,0].set_ylabel(y_label)
        axs[0,0].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        #------------------------------------------------------------------
    elif ii == 2:
        mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        pname = "Spring (MAM)"
        #axs[0,1].title.set_text("Spring (MAM)")
        #axs[0,1].title.set_text("Spring (Mar,Apr,May)")
        axs[0,1].set_ylabel(y_label)
        axs[0,1].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        #------------------------------------------------------------------
    elif ii == 3:
        mesh1 = axs[1,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        pname = "Summer (JJA)"
        #axs[1,0].title.set_text("Summer (JJA)")
        #axs[1,0].title.set_text("Summer (Jun,Jul,Aug)")
        axs[1,0].set_ylabel(y_label)
        axs[1,0].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        #------------------------------------------------------------------
    else:
        mesh1 = axs[1,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        pname = "Fall (SON)"
        #axs[1,1].title.set_text("Fall (SON)")
        #axs[1,1].title.set_text("Fall (Sep,Oct,Nov)")
        axs[1,1].set_ylabel(y_label)
        axs[1,1].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        #------------------------------------------------------------------

#stats_10 = list(range(int(np.shape(stats_vectors)[0])))
#stats_90 = list(range(int(np.shape(stats_vectors)[0])))

    #plt.colorbar(mesh1)

#fig.colorbar(mesh1, ax=axs.ravel().tolist())

cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist())
cbar.ax.set_ylabel(cbar_label, fontsize = cbar_fontSize)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.locator_params(nbins=cbar_nBins)

cbar_tick_labels_pre = [float(t.get_text().replace('−','-')) for t in cbar.ax.get_yticklabels()]
#cbar_tick_labels_pre = [float(t.get_text().replace('−','-')) for t in cbar.ax.get_yticklabels()]
#cbar_tick_labels_pre = [float(t.get_text()) for t in cbar.ax.get_yticklabels()]
#cbar_tick_labels_pre = [float(t.get_text() for t in cbar.ax.get_yticklabels())]

cbar_round = 4

cbar_tick_labels = [round(10**t, cbar_round) for t in cbar_tick_labels_pre]
#cbar_tick_labels = [10**t for t in cbar_tick_labels_pre]
cbar.set_ticklabels(cbar_tick_labels)


#fig.suptitle(f"$\it{fig_param_title}$ ~ {fig_main_title}")
#fig.suptitle(f"{fig_supTitle}", style = "italic")
#fig.suptitle(fig_supTitle)
fig.suptitle(fig_fullTitle)

#for a in fig.axes:
#    # Shrink the axes
#    box = a.get_position()
#    a.set_position([box.x0, box.y0, box.width * 0.9, box.height * 0.95])

#fig.text(s=f"{fig_param_title}", style = "italic",x=0.5, y=1.00, ha='center',va='center')
#fig.text(s=f"{fig_param_title}", style = "italic",x=0.5, y=.95, ha='center',va='center')
#fig.text(s=fig_main_title ,x=0.5, y=.90, ha='center',va='center')

plt.show()



