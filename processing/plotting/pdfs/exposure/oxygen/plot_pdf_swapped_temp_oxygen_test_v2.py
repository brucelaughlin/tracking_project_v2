# V2: Try the line plots for the day 0 data again

#------------------------------------------------------------------
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_oneFileTest_swapped.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_tenFileTest_swapped.p"
#pdf_file_name_pre = "pdf_data_output_releaseLoc_vs_settleTime_test3_swapped.p"
pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_run_test3_swapped.p"
#------------------------------------------------------------------

# -------------------------------------------------------------------
# Testing - choose 0-3 for test dex, corresponding to entries from 
#oxygen_limit_list = [2.2,3.1,4.1,6] # REMOVE THIS AFTER RUNNING CATCH SCRIPT AGAIN
test_dex = 2
#------------------------------------------------------------------


#pdf_file_name = pdf_file_name_pre[0:-2] + "_swapped.p"
pdf_file_name = pdf_file_name_pre
#------------------------------------------------------------------

# Titles, labels, random params
#------------------------------------------------------------------
fig_param_title = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD\n"
#fig_param_title = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"
#fig_main_title_pre = "PDFs of cumulative number of days larvae were exposed to DO2 concentrations below {} mg/L (y-axis) vs Release Location (x-axis).\nGrouped according to season of release."






y_label_pre = "exposure time (days)"
#y_label_pre = "exposure time to [DO2] below {} mg/L (days)"

day_0_title_pre = "probability of 0 days of exposure"
#day_0_title_pre = "probability of 0 days of exposure to [DO2] below {} mg/L"

day_1_60_title_pre = "probability of 1-60 days of exposure"
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


file = open(pdf_raw_file,'rb')

pdf_list_exposure_T_source_swapped,pdf_list_of_lists_O2_source_swapped,pdf_list_connectivity_swapped,pdf_list_settleTime_swapped,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels,first_continent_box_dex,oxygen_limit_list = pickle.load(file)
#pdf_list_exposure_T_source_swapped,pdf_list_of_lists_O2_source_swapped,pdf_list_connectivity_swapped,pdf_list_settleTime_swapped,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels,first_continent_box_dex = pickle.load(file)
file.close()


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Testing - choose 0-3 for test dex, corresponding to entries from 
#oxygen_limit_list = [2.2,3.1,4.1,6] # REMOVE THIS AFTER RUNNING CATCH SCRIPT AGAIN
#test_dex = 3
# -------------------------------------------------------------------
pdf_list_exposure_oxygen_source = pdf_list_of_lists_O2_source_swapped[test_dex]
# -------------------------------------------------------------------
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Dynamic title and labels
#fig_main_title = fig_main_title_pre.format(oxygen_limit_list[test_dex])

#fig_main_title = fig_main_title_pre.format("$\bf{{{x}}}$".format(x=oxygen_limit_list[test_dex]))
#fig_main_title = "$\bf{{{x}}}$".format(x=oxygen_limit_list[test_dex])

fig_main_title_pre_p1 = "PDFs of cumulative number of days settled larvae were exposed to DO2 concentrations below "
fig_main_title_pre_p2 = r'$\bf{{{a}}}$'.format(a=oxygen_limit_list[test_dex])
fig_main_title_pre_p3 = " mg/L (y-axis) vs Release Location (x-axis).\nGrouped according to season of release.\n(Day 1 PDF plotted below main plot as a line plot)"
#fig_main_title_pre_p3 = " mg/L (y-axis) vs Release Location (x-axis).\nGrouped according to season of release."

fig_main_title = fig_main_title_pre_p1 + fig_main_title_pre_p2 + fig_main_title_pre_p3



#plt.rc('text', usetex=True)
#fig_main_title = r"This is text. I want the word '\textbf{blah}' to be Bold."



y_label = y_label_pre.format(oxygen_limit_list[test_dex])


# -------------------------------------------------------------------
fig_fullTitle = fig_main_title + "\n" + fig_param_title
# -------------------------------------------------------------------




pdf_list_full = []
pdf_list_line = []


# MESH DATA
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
    pdf_full = pdf_full[:,1:]
    pdf_list_full.append(pdf_full)
    if np.amax(pdf_full) > pdf_max_val:
        pdf_max_val = np.amax(pdf_full)
    if np.amin(np.ma.masked_invalid(pdf_full)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf_full))
# -------------------------------------------------------------------

# LINE DATA
pdf_max_val_day1 = -999999
pdf_min_val_day1 = 999999
# -------------------------------------------------------------------
for pdf in pdf_list_exposure_oxygen_source[1:]:
    pdf_full = pdf.copy()
    # -------------------------------------------------------------------
    # Convert to probability
    row_sums = pdf_full.sum(axis=1)
    pdf_full = pdf_full / row_sums[:, np.newaxis]
    #pdf_full = np.log10(pdf_full)
    # -------------------------------------------------------------------
    pdf_full = pdf_full[:,0]
    pdf_list_line.append(pdf_full)
    if np.amax(pdf_full) > pdf_max_val_day1:
        pdf_max_val_day1 = np.amax(pdf_full)
    if np.amin(np.ma.masked_invalid(pdf_full)) < pdf_min_val_day1:
        pdf_min_val_day1 = np.amin(np.ma.masked_invalid(pdf_full))
# -------------------------------------------------------------------





# Determined elsewhere (see/run "check_box_numbers.py")
#first_continent_box_dex = 17
#first_continent_box_dex = 20
num_dummy_lines = 1

n_rows = int(np.shape(pdf_list_exposure_oxygen_source[0])[0]) + num_dummy_lines
n_columns = int(np.shape(pdf_list_exposure_oxygen_source[0])[1])
X = np.arange(-0.5, n_rows, 1)
Y = np.arange(0.5, n_columns, 1)
#Y = np.arange(-0.5, n_columns, 1)


#y_tick_day1_spacing = .15
y_tick_day1_number = 4
y_ticks_day1 = list(np.linspace(pdf_min_val_day1,pdf_max_val_day1, y_tick_day1_number))
y_tick_day1_round = 3
y_ticks_day1 = [round(t,y_tick_day1_round) for t in y_ticks_day1]


#cbar_tick_labels_pre = [float(t.get_text().replace('−','-')) for t in cbar.ax.get_yticklabels()]
#cbar_round = 4
#cbar_tick_labels = [round(10**t, cbar_round) for t in cbar_tick_labels_pre]
#cbar.set_ticklabels(cbar_tick_labels)

day_1_twin_label = "0 days"
days_fullM1_twin_label = "1-60 days"

day_0_title = day_0_title_pre.format(oxygen_limit_list[test_dex])
day_1_60_title_pre2 = day_1_60_title_pre.format(oxygen_limit_list[test_dex])

seasons = ["(Dec, Jan, Feb)","(Mar, Apr, May)","(Jun, Jul, Aug)","(Sep, Oct, Nov)"]
day_1_60_titles = []
for ii in range(len(seasons)):
    day_1_60_titles.append("{}     {}".format(day_1_60_title_pre2,seasons[ii]))

#extra_ticks = [1]
#axs[0,0].set_yticks(list(plt.xticks()[0]) + extra_ticks)


v_scale = 6

fig,axs = plt.subplots(4,2, height_ratios = [v_scale,1,v_scale,1], constrained_layout=True)
#fig,axs = plt.subplots(2,2)

label_fontSize = 10
#label_fontSize = 8

for ii in range(len(pdf_list_full)):

    # 2D full plot (minus day 1)
    pdf_plot = pdf_list_full[ii]
    pdf_separated = np.empty((np.shape(pdf_plot)[0] + num_dummy_lines,np.shape(pdf_plot)[1]))
    pdf_separated[:] = np.nan
    pdf_separated[0:first_continent_box_dex,:] = pdf_plot[0:first_continent_box_dex,:]
    pdf_separated[first_continent_box_dex + num_dummy_lines:,:] = pdf_plot[first_continent_box_dex:,:]


    # 1D day one pdf
    pdf_day1 = pdf_list_line[ii]
    pdf_day1_separated = np.empty((np.shape(pdf_day1)[0] + num_dummy_lines))
    pdf_day1_separated[:] = np.nan
    pdf_day1_separated[0:first_continent_box_dex] = pdf_day1[0:first_continent_box_dex]
    pdf_day1_separated[first_continent_box_dex + num_dummy_lines:] = pdf_day1[first_continent_box_dex:]





    if ii == 0:
    #if ii == 1:
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,0].title.set_text(day_1_60_titles[ii])
        #axs[0,0].title.set_text("Winter (Dec,Jan,Feb)")
        axs[0,0].set_ylabel(y_label)
        axs[0,0].yaxis.label.set(fontsize=15)
        axs[1,0].plot(pdf_day1_separated)
        axs[1,0].set_ylim([pdf_min_val_day1,pdf_max_val_day1])
        axs[1,0].margins(x=0)
        axs[1,0].set_yticks(y_ticks_day1)
        axs[1,0].yaxis.grid(True)
        axs[1,0].title.set_text(day_0_title)
        axs[1,0].set_ylabel(cbar_label)
        #------------------------------------------------------------------
        #axs[0,0].set_xticks(tick_positions)
        #axs[0,0].set_xticklabels(tick_labels, fontsize=label_fontSize)
        ##axs[0,0].set_xticklabels(tick_labels_double_X, fontsize=label_fontSize)
        #axs[1,0].set_xticks(tick_positions)
        #axs[1,0].set_xticklabels(tick_labels_single_X, fontsize=label_fontSize)
        #------------------------------------------------------------------
    elif ii == 1:
    #elif ii == 2:
        mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,1].title.set_text(day_1_60_titles[ii])
        #axs[0,1].title.set_text("Spring (Mar,Apr,May)")
        axs[0,1].set_ylabel(y_label)
        axs[0,1].yaxis.label.set(fontsize=15)
        axs[1,1].plot(pdf_day1_separated)
        axs[1,1].set_ylim([pdf_min_val_day1,pdf_max_val_day1])
        axs[1,1].margins(x=0)
        axs[1,1].set_yticks(y_ticks_day1)
        axs[1,1].yaxis.grid(True)
        axs[1,1].title.set_text(day_0_title)
        axs[1,1].set_ylabel(cbar_label)
        #------------------------------------------------------------------
        #ax2 = axs[0,1].twinx()
        #ax2.set_ylabel(days_fullM1_twin_label, fontsize=15)
        #ax2 = axs[1,1].twinx()
        #ax2.set_ylabel(day_1_twin_label, fontsize=15)
        #------------------------------------------------------------------
        #axs[0,1].set_xticks(tick_positions)
        #axs[0,1].set_xticklabels(tick_labels, fontsize=label_fontSize)
        ##axs[0,1].set_xticklabels(tick_labels_double_X, fontsize=label_fontSize)
        #axs[1,1].set_xticks(tick_positions)
        #axs[1,1].set_xticklabels(tick_labels_single_X, fontsize=label_fontSize)
        #------------------------------------------------------------------
    elif ii == 2:
    #elif ii == 3:
        mesh1 = axs[2,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[2,0].title.set_text(day_1_60_titles[ii])
        #axs[2,0].title.set_text("Summer (Jun,Jul,Aug)")
        axs[2,0].set_ylabel(y_label)
        axs[2,0].yaxis.label.set(fontsize=15)
        axs[3,0].plot(pdf_day1_separated)
        axs[3,0].set_ylim([pdf_min_val_day1,pdf_max_val_day1])
        axs[3,0].margins(x=0)
        axs[3,0].set_yticks(y_ticks_day1)
        axs[3,0].yaxis.grid(True)
        axs[3,0].title.set_text(day_0_title)
        axs[3,0].set_ylabel(cbar_label)
        #------------------------------------------------------------------
        #axs[2,0].set_xticks(tick_positions)
        #axs[2,0].set_xticklabels(tick_labels, fontsize=label_fontSize)
        ##axs[2,0].set_xticklabels(tick_labels_double_X, fontsize=label_fontSize)
        #axs[3,0].set_xticks(tick_positions)
        #axs[3,0].set_xticklabels(tick_labels_single_X, fontsize=label_fontSize)
        #------------------------------------------------------------------
    else:
        mesh1 = axs[2,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[2,1].title.set_text(day_1_60_titles[ii])
        #axs[2,1].title.set_text("Fall (Sep,Oct,Nov)")
        axs[2,1].set_ylabel(y_label)
        axs[2,1].yaxis.label.set(fontsize=15)
        axs[3,1].plot(pdf_day1_separated)
        axs[3,1].set_ylim([pdf_min_val_day1,pdf_max_val_day1])
        axs[3,1].margins(x=0)
        axs[3,1].set_yticks(y_ticks_day1)
        axs[3,1].yaxis.grid(True)
        axs[3,1].title.set_text(day_0_title)
        axs[3,1].set_ylabel(cbar_label)
        #------------------------------------------------------------------
        #ax2 = axs[2,1].twinx()
        #ax2.set_ylabel(days_fullM1_twin_label, fontsize=15)
        #ax2 = axs[3,1].twinx()
        #ax2.set_ylabel(day_1_twin_label, fontsize=15)
        #------------------------------------------------------------------
        #axs[2,1].set_xticks(tick_positions)
        #axs[2,1].set_xticklabels(tick_labels, fontsize=label_fontSize)
        ##axs[2,1].set_xticklabels(tick_labels_double_X, fontsize=label_fontSize)
        #axs[3,1].set_xticks(tick_positions)
        #axs[3,1].set_xticklabels(tick_labels_single_X, fontsize=label_fontSize)
        #------------------------------------------------------------------





#    if ii == 1:
#        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
#        #axs[0,0].title.set_text("Winter (Dec,Jan,Feb)")
#        axs[0,0].set_ylabel(y_label)
#        axs[0,0].yaxis.label.set(fontsize=15)
#        #------------------------------------------------------------------
#        #------------------------------------------------------------------
#    elif ii == 2:
#        mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
#        #axs[0,1].title.set_text("Spring (Mar,Apr,May)")
#        axs[0,1].set_ylabel(y_label)
#        axs[0,1].yaxis.label.set(fontsize=15)
#        #------------------------------------------------------------------
#        #------------------------------------------------------------------
#    elif ii == 3:
#        mesh1 = axs[1,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
#        #axs[1,0].title.set_text("Summer (Jun,Jul,Aug)")
#        axs[1,0].set_ylabel(y_label)
#        axs[1,0].yaxis.label.set(fontsize=15)
#        #------------------------------------------------------------------
#        #------------------------------------------------------------------
#    else:
#        mesh1 = axs[1,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
#        #axs[1,1].title.set_text("Fall (Sep,Oct,Nov)")
#        axs[1,1].set_ylabel(y_label)
#        axs[1,1].yaxis.label.set(fontsize=15)
#        #------------------------------------------------------------------
#        #------------------------------------------------------------------






#------------------------------------------------------------------
# Colorbar
#------------------------------------------------------------------
cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist())
cbar.ax.set_ylabel(cbar_label, fontsize = cbar_fontSize)
cbar.ax.yaxis.set_label_position('left')
cbar.ax.locator_params(nbins=cbar_nBins)
cbar_tick_labels_pre = [float(t.get_text().replace('−','-')) for t in cbar.ax.get_yticklabels()]
cbar_round = 4
cbar_tick_labels = [round(10**t, cbar_round) for t in cbar_tick_labels_pre]
cbar.set_ticklabels(cbar_tick_labels)
#------------------------------------------------------------------


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



