#V4: just do simple colorbar, raw probablity values

#------------------------------------------------------------------
#pdf_file_name_pre = "pdf_data_output_seasonal_v3_tenFileTest.p"
pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_run_test3_swapped.p"
#------------------------------------------------------------------

#------------------------------------------------------------------
#pdf_file_name = pdf_file_name_pre[0:-2] + "_swapped.p"
pdf_file_name = pdf_file_name_pre
#------------------------------------------------------------------


#------------------------------------------------------------------
fig_paramTitle = "wc_15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"
fig_mainTitle = "PDFs of average temperature experienced (y-axis) vs release box number (x-axis).\nGrouped according to season of release."
#fig_mainTitle = "PDFs of Average Temperature Experienced (y-axis) vs Release Box Number (x-axis).\nGrouped according to season of release."
#fig_mainTitle = "PDFs of Average Temperature Experienced (y-axis) vs Release Location (x-axis).\nGrouped according to season of release."

fig_fullTitle = fig_mainTitle + "\n" + fig_paramTitle

yLabel = "average temperature experienced (C$^\circ$)"



cbar_label = "probability"
#cbar_label = "Log base 10 of probability"
#cbar_label = "Log base 10 value of probability"
cbar_fontSize = 15
cbar_nBins = 20
#cbar_nBins = 25

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


# Trim histograms... done by eye by plotting this same plot...
T_min = 5
T_max = 25
n_ticks = int((T_max-T_min)/2) - 1 #kinda a hack, seems to match the naive plot without specified ticks

bin_size = 0.1
bin_T_min = int(T_min/bin_size)
bin_T_max = int(T_max/bin_size)
bin_T_range = bin_T_max - bin_T_min

#y_ticks = list(range(0,bin_T_range+1,y_tick_step))
#y_ticks = np.linspace(bin_T_min,bin_T_max,n_ticks)
y_ticks = np.linspace(0,bin_T_range,n_ticks)
y_ticks = [int(x) for x in y_ticks]

y_tick_labels =  np.linspace(T_min,T_max,n_ticks)
y_tick_labels = [str(x) for x in y_tick_labels]



file = open(pdf_raw_file,'rb')
pdf_list_exposure_T_source_swapped,pdf_list_of_lists_O2_source_swapped,pdf_list_connectivity_swapped,pdf_list_settleTime_swapped,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels,first_continent_box_dex,oxygen_limit_list = pickle.load(file)
#pdf_list_exposure_T_source,pdf_list_exposure_oxygen_source,pdf_list_connectivity,pdf_list_settleTime,settlement_boxes_test_array,settlement_times_test_array,counter_array = pickle.load(file)
#pdf_list_exposure_T_source,pdf_list_exposure_oxygen_source,pdf_list_connectivity,pdf_list_settleTime,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels = pickle.load(file)
file.close()



pdf_list_full = []

pdf_max_val = -999999
pdf_min_val = 999999

# -------------------------------------------------------------------
for pdf in pdf_list_exposure_T_source_swapped[1:]:
#for pdf in pdf_list_exposure_T_source[1:]:
    pdf_full = pdf.copy()

    pdf_full = pdf_full[:,bin_T_min:bin_T_max+1]


    row_sums = pdf_full.sum(axis=1)
    pdf_full = pdf_full / row_sums[:, np.newaxis]
    pdf_full = np.log10(pdf_full)
    pdf_list_full.append(pdf_full)
    if np.amax(pdf_full) > pdf_max_val:
        pdf_max_val = np.amax(pdf_full)
    if np.amin(np.ma.masked_invalid(pdf_full)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf_full))
# -------------------------------------------------------------------


num_dummy_lines = 1

n_rows = int(np.shape(pdf_list_full[0])[0]) + num_dummy_lines
n_columns = int(np.shape(pdf_list_full[0])[1])
#n_rows = int(np.shape(pdf_list_exposure_T_source_swapped[0])[0]) + num_dummy_lines
#n_columns = int(np.shape(pdf_list_exposure_T_source_swapped[0])[1])
X = np.arange(-0.5, n_rows, 1)
Y = np.arange(-0.5, n_columns, 1)


#Modify "tick_positions" for plotting (we add an empty bin before the first continent box)
tick_positions_modified = tick_positions.copy()
for ii in range(len(tick_positions_modified)):
    if tick_positions_modified[ii] >= first_continent_box_dex:
        tick_positions_modified[ii] += num_dummy_lines


stagger_dex = 0
tick_labels_double_X = []
for ii in range(len(tick_labels)):
    stagger_dex += 1
    if (tick_positions[ii] < first_continent_box_dex) and (stagger_dex % 2 == 0):
    #if (tick_positions[ii]+1 >= 11) and (tick_positions[ii]+1 <= 17) and (stagger_dex % 2 == 0):
        #tick_labels_double_X.append("{}\n\n{}".format(tick_positions[ii],tick_labels[ii]))
        tick_labels_double_X.append("{}\n\n{}".format(tick_positions[ii]+1,tick_labels[ii]))
    else:
        #tick_labels_double_X.append("{}\n{}".format(tick_positions[ii],tick_labels[ii]))
        tick_labels_double_X.append("{}\n{}".format(tick_positions[ii]+1,tick_labels[ii]))




seasons = ["(Dec, Jan, Feb)","(Mar, Apr, May)","(Jun, Jul, Aug)","(Sep, Oct, Nov)"]
#season_title_pre = ""
season_titles = []
for ii in range(len(seasons)):
    season_titles.append("{}".format(seasons[ii]))


#x_line_positions = []
#for tick in tick_positions:
#    if tick < first_continent_box_dex:
#        x_line_positions.append(tick-0.5)
#    x_line_positions.append(tick-0.5)
    
x_line_positions = []
for ii in range(np.shape(pdf_list_full[0])[0]):
#for ii in range(2,np.shape(pdf_list_full[0])[0],2):
    #if ii < first_continent_box_dex:
    if ii < first_continent_box_dex:
        if ii%2 == 0:
            x_line_positions.append(ii+1.5)
            #x_line_positions.append(ii+0.5)
    else:
        x_line_positions.append(ii+0.5)
        #x_line_positions.append(ii+1.5)



v_scale = 6

fig,axs = plt.subplots(2,2)


plt.setp(axs,xticks=tick_positions_modified,xticklabels=tick_labels_double_X)
#plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels_double_X)


label_fontSize = 10
#label_fontSize = 8

for ii in range(len(pdf_list_full)):
#for ii in range(1,len(pdf_list_exposure_T_source)):

    #pdf_plot_pre1 = pdf_list_exposure_T_source[ii]
    
    #row_sums = pdf_plot_pre1.sum(axis=1)
    #pdf_plot_pre2 = pdf_plot_pre1 / row_sums[:, np.newaxis]

    # 2D full plot (minus day 1)
    pdf_plot = pdf_list_full[ii]
    pdf_separated = np.empty((np.shape(pdf_plot)[0] + num_dummy_lines,np.shape(pdf_plot)[1]))
    pdf_separated[:] = np.nan
    pdf_separated[0:first_continent_box_dex,:] = pdf_plot[0:first_continent_box_dex,:]
    pdf_separated[first_continent_box_dex + num_dummy_lines:,:] = pdf_plot[first_continent_box_dex:,:]
    #pdf_separated[0:first_continent_box_dex-1,:] = pdf_plot[0:first_continent_box_dex-1,:]
    #pdf_separated[first_continent_box_dex-1 + num_dummy_lines:,:] = pdf_plot[first_continent_box_dex-1:,:]



    if ii == 0:
    #if ii == 1:
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,0].title.set_text(season_titles[ii])
        axs[0,0].set_yticks(y_ticks)
        axs[0,0].set_yticklabels(y_tick_labels)
        axs[0,0].set_ylabel(yLabel)
        axs[0,0].yaxis.label.set(fontsize=15)
        axs[0,0].set_xticks(x_line_positions,minor=True)
        axs[0,0].xaxis.grid(True,which='minor')
        #------------------------------------------------------------------
#        axs[0,0].set_xticks(tick_positions)
#        axs[0,0].set_xticklabels(tick_labels_double_X, fontsize=label_fontSize)
        #------------------------------------------------------------------
    elif ii == 1:
    #elif ii == 2:
        mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,1].title.set_text(season_titles[ii])
        axs[0,1].set_yticks(y_ticks)
        axs[0,1].set_yticklabels(y_tick_labels)
        axs[0,1].set_xticks(x_line_positions,minor=True)
        axs[0,1].xaxis.grid(True,which='minor')
        #axs[0,1].set_ylabel(yLabel)
        #axs[0,1].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        #------------------------------------------------------------------
    elif ii == 2:
    #elif ii == 3:
        mesh1 = axs[1,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,0].title.set_text(season_titles[ii])
        axs[1,0].set_yticks(y_ticks)
        axs[1,0].set_yticklabels(y_tick_labels)
        axs[1,0].set_ylabel(yLabel)
        axs[1,0].yaxis.label.set(fontsize=15)
        axs[1,0].set_xticks(x_line_positions,minor=True)
        axs[1,0].xaxis.grid(True,which='minor')
        #------------------------------------------------------------------
        #------------------------------------------------------------------
    else:
        mesh1 = axs[1,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,1].title.set_text(season_titles[ii])
        axs[1,1].set_yticks(y_ticks)
        axs[1,1].set_yticklabels(y_tick_labels)
        axs[1,1].set_xticks(x_line_positions,minor=True)
        axs[1,1].xaxis.grid(True,which='minor')
        #axs[1,1].set_ylabel(yLabel)
        #axs[1,1].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        #------------------------------------------------------------------


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




#fig.suptitle(f"$\it{fig_paramTitle}$ ~ {fig_mainTitle}")
#fig.suptitle(f"{fig_supTitle}", style = "italic")
#fig.suptitle(fig_supTitle)
fig.suptitle(fig_fullTitle)

#for a in fig.axes:
#    # Shrink the axes
#    box = a.get_position()
#    a.set_position([box.x0, box.y0, box.width * 0.9, box.height * 0.95])

#fig.text(s=f"{fig_paramTitle}", style = "italic",x=0.5, y=1.00, ha='center',va='center')
#fig.text(s=f"{fig_paramTitle}", style = "italic",x=0.5, y=.95, ha='center',va='center')
#fig.text(s=fig_mainTitle ,x=0.5, y=.90, ha='center',va='center')

plt.show()



