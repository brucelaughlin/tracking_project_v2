# v1: Now try presenting as cdf


pld_day_first = 31
#pld_day_first = 61
#pld_day_first = 91

# Define run length (ie last day of PLD)
#pld_day_last = 151
#pld_day_last = 31
#pld_day_last = 61
pld_day_last = 91


y_ticks_full_spacing = 10


# The extra one is an artifact from the output (we get data at time 0)
pld_day_first = pld_day_first - 1
pld_day_last = pld_day_last - 1



#------------------------------------------------------------------
pdf_file_name_pre = "pdf_data_output_releaseLoc_vs_settleTime_test3.p"
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_test4_physics_only_AKs_1en5.p"
#------------------------------------------------------------------


#------------------------------------------------------------------
pdf_file_name = pdf_file_name_pre[0:-2] + "_swapped.p"
#------------------------------------------------------------------


#------------------------------------------------------------------
fig_paramTitle = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, " + "PLD days {}-{}".format(pld_day_first,pld_day_last)
#fig_paramTitle = "wc_15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"
fig_mainTitle = "CDFs of Settlement Day within PLD (y-axis) vs Release Location (x-axis).\nGrouped according to season of release.\n(Day 1 probability plotted below main plot as a line plot)"

fig_fullTitle = fig_mainTitle + "\n" + fig_paramTitle

#------------------------------------------------------------------


#------------------------------------------------------------------
# Plots are dominated by settlement on day 1 of settlement window - so here's a switch to remove day 1 from the plots
#switch_day_one = True
switch_day_one = False
#------------------------------------------------------------------


import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import datetime
import cm_ocean

yLabel = "Days {}-{}".format(pld_day_first+1,pld_day_last)
#yLabel = "Days since end of PLD"
cbar_label = "probability"


day_0_title = "probability of settlement on day {}".format(pld_day_first)
day_1_60_title_pre = "CDFs for settlement on days {}-{}  (settler release months:".format(pld_day_first+1,pld_day_last)
#day_1_60_title_pre = "probability of settlement on day 2-  (settler release months:"



#------------------------------------------------------------------
base_path = '/home/blaughli/tracking_project/'
pdf_raw_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/'

#pdf_raw_file = pdf_raw_directory + 'x_drift_30_noSwim/' + pdf_file_name
pdf_raw_file = pdf_raw_directory + pdf_file_name
#------------------------------------------------------------------

#base_datetime = datetime.datetime(1970,1,1,0,0,0)


file = open(pdf_raw_file,'rb')
pdf_list_connectivity,pdf_list_settleTime,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels = pickle.load(file)
file.close()

#pdf_conn_list = []
#for pdf in pdf_list_connectivity:
#    pdf_conn_list.append(np.log10(np.transpose(pdf)))


cdf_list_full = []
cdf_list_line = []

cdf_max_val = -999999
cdf_min_val = 999999
cdf_max_val_day1 = -999999
cdf_min_val_day1 = 999999

for pdf in pdf_list_settleTime[1:]:
    pdf_time_full = pdf
    row_sums = pdf_time_full.sum(axis=1)
    pdf_time_full = pdf_time_full / row_sums[:, np.newaxis]
#    pdf_time_full = np.log10(pdf_time_full)

    #cdf_time_full = pdf_time_full.copy()
    cdf_time_full = np.zeros_like(pdf_time_full)
    for ii in range(np.shape(cdf_time_full)[1]):
        cdf_time_full[:,ii] = np.sum(pdf_time_full[:,0:ii+1],axis=1)



    cdf_time_line = cdf_time_full[:,0]
    cdf_list_line.append(cdf_time_line)
    
    cdf_time_full = cdf_time_full[:,1:]
    cdf_list_full.append(cdf_time_full)

    if np.amax(cdf_time_full) > cdf_max_val:
        cdf_max_val = np.amax(cdf_time_full)
    if np.amin(np.ma.masked_invalid(cdf_time_full)) < cdf_min_val:
        cdf_min_val = np.amin(np.ma.masked_invalid(cdf_time_full))
    
    if np.amax(cdf_time_line) > cdf_max_val_day1:
        cdf_max_val_day1 = np.amax(cdf_time_line)
    if np.amin(np.ma.masked_invalid(cdf_time_line)) < cdf_min_val_day1:
        cdf_min_val_day1 = np.amin(np.ma.masked_invalid(cdf_time_line))



## Also need min/max for day1 pdfs
#pdf_max_val_day1 = -999999
#pdf_min_val_day1 = 999999
#for pdf in pdf_list_settleTime[1:]:
#    pdf_time_full = pdf
#    row_sums = pdf_time_full.sum(axis=1)
#    pdf_time_full = pdf_time_full / row_sums[:, np.newaxis]
##    pdf_time_full = np.log10(pdf_time_full)
#    pdf_time_full = pdf_time_full[:,0]
#    pdf_list_line.append(pdf_time_full)
#    if np.amax(pdf_time_full) > pdf_max_val_day1:
#        pdf_max_val_day1 = np.amax(pdf_time_full)
#    if np.amin(np.ma.masked_invalid(pdf_time_full)) < pdf_min_val_day1:
#        pdf_min_val_day1 = np.amin(np.ma.masked_invalid(pdf_time_full))



# Determined elsewhere (see/run "check_box_numbers.py")
first_continent_box_dex = 20
num_dummy_lines = 1


n_boxes_seeded = int(np.shape(pdf_list_settleTime[0])[0]) + num_dummy_lines
n_times = int(np.shape(pdf_list_settleTime[0])[1])
X = np.arange(-0.5, n_boxes_seeded, 1)
if switch_day_one:
    Y = np.arange(-0.5, n_times, 1)
else:
    Y = np.arange(0.5, n_times, 1)

v_scale = 6



stagger_dex = 0
tick_labels_double_X = []
for ii in range(len(tick_labels)):
    stagger_dex += 1
    if (tick_positions[ii]+1 >= 11) and (tick_positions[ii]+1 <= 17) and (stagger_dex % 2 == 0):
        tick_labels_double_X.append("{}\n\n{}".format(tick_positions[ii]+1,tick_labels[ii]))
    else:
        tick_labels_double_X.append("{}\n{}".format(tick_positions[ii]+1,tick_labels[ii]))

tick_labels_single_X = [x + 1 for x in tick_positions]

y_tick_day1_spacing = .1
#y_tick_day1_spacing = .15
y_ticks_day1 = np.arange(cdf_min_val_day1,cdf_max_val_day1, y_tick_day1_spacing)
y_tick_day1_round = 2
y_ticks_day1 = [round(t,y_tick_day1_round) for t in y_ticks_day1]



#y_ticks_full_spacing = int((pld_day_last-pld_day_first)/5) 
#y_tick_labels_full = np.arange(pld_day_first-1,pld_day_last-1,y_ticks_full_spacing)
#y_tick_labels_full = np.arange(pld_day_first,pld_day_last,y_ticks_full_spacing)
y_tick_labels_full = np.arange(pld_day_first,pld_day_last,y_ticks_full_spacing) + 1
y_ticks_full = y_tick_labels_full - pld_day_first
#y_ticks_full = y_tick_labels_full - pld_day_first + 1

print(y_tick_labels_full)
print(y_ticks_full)

seasons = ["Dec, Jan, Feb)","Mar, Apr, May)","Jun, Jul, Aug)","Sep, Oct, Nov)"]
day_1_60_titles = []
for ii in range(len(seasons)):
    day_1_60_titles.append("{} {}".format(day_1_60_title_pre,seasons[ii]))
    #day_1_60_titles.append("{} {}".format(day_1_60_title_pre2,seasons[ii]))





fig,axs = plt.subplots(4,2, height_ratios = [v_scale,1,v_scale,1], constrained_layout=True)


label_fontSize = 10
#label_fontSize = 8

for ii in range(len(cdf_list_full)):
#for ii in range(1,len(cdf_list_settleTime)):

    #cdf_plot_pre1 = cdf_list_settleTime[ii]
    
    #row_sums = cdf_plot_pre1.sum(axis=1)
    #cdf_plot_pre2 = cdf_plot_pre1 / row_sums[:, np.newaxis]

    # 2D full plot (minus day 1)
    cdf_plot = cdf_list_full[ii]
    #cdf_plot = cdf_plot_pre2[:,1:]
    cdf_separated = np.empty((np.shape(cdf_plot)[0] + num_dummy_lines,np.shape(cdf_plot)[1]))
    cdf_separated[:] = np.nan
    cdf_separated[0:first_continent_box_dex,:] = cdf_plot[0:first_continent_box_dex,:]
    cdf_separated[first_continent_box_dex + num_dummy_lines:,:] = cdf_plot[first_continent_box_dex:,:]

    # 1D day one cdf
    cdf_day1 = cdf_list_line[ii]
    #cdf_day1 = cdf_plot_pre2[:,0]
    cdf_day1_separated = np.empty((np.shape(cdf_day1)[0] + num_dummy_lines))
    cdf_day1_separated[:] = np.nan
    cdf_day1_separated[0:first_continent_box_dex] = cdf_day1[0:first_continent_box_dex]
    cdf_day1_separated[first_continent_box_dex + num_dummy_lines:] = cdf_day1[first_continent_box_dex:]

    if ii == 0:
    #if ii == 1:
        mesh1 = axs[0,0].pcolormesh(X,Y,cdf_separated.T,cmap='jet',vmin=cdf_min_val,vmax=cdf_max_val)
        axs[0,0].title.set_text(day_1_60_titles[ii])
        axs[0,0].set_ylabel(yLabel)
        axs[0,0].yaxis.label.set(fontsize=15)
        axs[0,0].set_yticks(y_ticks_full)
        axs[0,0].set_yticklabels(y_tick_labels_full)
        axs[0,0].set_xticks(tick_positions)
        axs[0,0].set_xticklabels(tick_labels, fontsize=label_fontSize)
        #------------------------------------------------------------------
        axs[1,0].plot(cdf_day1_separated)
        axs[1,0].set_ylim([cdf_min_val_day1,cdf_max_val_day1])
        axs[1,0].margins(x=0)
        axs[1,0].set_yticks(y_ticks_day1)
        axs[1,0].set_ylabel(cbar_label)
        axs[1,0].yaxis.grid(True)
        axs[1,0].set_xticks(tick_positions)
        axs[1,0].set_xticklabels(tick_labels_single_X, fontsize=label_fontSize)
        axs[1,0].title.set_text(day_0_title)
        #------------------------------------------------------------------
    elif ii == 1:
    #elif ii == 2:
        mesh1 = axs[0,1].pcolormesh(X,Y,cdf_separated.T,cmap='jet',vmin=cdf_min_val,vmax=cdf_max_val)
        axs[0,1].title.set_text(day_1_60_titles[ii])
        axs[0,1].set_ylabel(yLabel)
        axs[0,1].yaxis.label.set(fontsize=15)
        axs[0,1].set_yticks(y_ticks_full)
        axs[0,1].set_yticklabels(y_tick_labels_full)
        axs[0,1].set_xticks(tick_positions)
        axs[0,1].set_xticklabels(tick_labels, fontsize=label_fontSize)
        #------------------------------------------------------------------
        axs[1,1].plot(cdf_day1_separated)
        axs[1,1].set_ylim([cdf_min_val_day1,cdf_max_val_day1])
        axs[1,1].margins(x=0)
        axs[1,1].set_yticks(y_ticks_day1)
        axs[1,1].set_ylabel(cbar_label)
        axs[1,1].yaxis.grid(True)
        axs[1,1].set_xticks(tick_positions)
        axs[1,1].set_xticklabels(tick_labels_single_X, fontsize=label_fontSize)
        axs[1,1].title.set_text(day_0_title)
        #------------------------------------------------------------------
    elif ii == 2:
    #elif ii == 3:
        mesh1 = axs[2,0].pcolormesh(X,Y,cdf_separated.T,cmap='jet',vmin=cdf_min_val,vmax=cdf_max_val)
        axs[2,0].set_ylabel(yLabel)
        axs[2,0].yaxis.label.set(fontsize=15)
        axs[2,0].set_yticks(y_ticks_full)
        axs[2,0].set_yticklabels(y_tick_labels_full)
        axs[2,0].set_xticks(tick_positions)
        axs[2,0].set_xticklabels(tick_labels, fontsize=label_fontSize)
        axs[2,0].title.set_text(day_1_60_titles[ii])
        #------------------------------------------------------------------
        axs[3,0].plot(cdf_day1_separated)
        axs[3,0].set_ylim([cdf_min_val_day1,cdf_max_val_day1])
        axs[3,0].margins(x=0)
        axs[3,0].set_yticks(y_ticks_day1)
        axs[3,0].set_ylabel(cbar_label)
        axs[3,0].yaxis.grid(True)
        axs[3,0].set_xticks(tick_positions)
        axs[3,0].set_xticklabels(tick_labels_single_X, fontsize=label_fontSize)
        axs[3,0].title.set_text(day_0_title)
        #------------------------------------------------------------------
    else:
        mesh1 = axs[2,1].pcolormesh(X,Y,cdf_separated.T,cmap='jet',vmin=cdf_min_val,vmax=cdf_max_val)
        axs[2,1].title.set_text(day_1_60_titles[ii])
        axs[2,1].set_ylabel(yLabel)
        axs[2,1].yaxis.label.set(fontsize=15)
        axs[2,1].set_yticks(y_ticks_full)
        axs[2,1].set_yticklabels(y_tick_labels_full)
        axs[2,1].set_xticks(tick_positions)
        axs[2,1].set_xticklabels(tick_labels, fontsize=label_fontSize)
        #------------------------------------------------------------------
        axs[3,1].plot(cdf_day1_separated)
        axs[3,1].set_ylim([cdf_min_val_day1,cdf_max_val_day1])
        axs[3,1].margins(x=0)
        axs[3,1].set_yticks(y_ticks_day1)
        axs[3,1].set_ylabel(cbar_label)
        axs[3,1].yaxis.grid(True)
        axs[3,1].set_xticks(tick_positions)
        axs[3,1].set_xticklabels(tick_labels_single_X, fontsize=label_fontSize)
        axs[3,1].title.set_text(day_0_title)
        #------------------------------------------------------------------

    #plt.colorbar(mesh1)

#fig.colorbar(mesh1, ax=axs.ravel().tolist())



cbar_label = "cumulative probability"
cbar_fontSize = 15
cbar_nBins = 20

#cbaxes = fig.add_axes([0.95,0.1, 0.03, 0.8])

#cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist(), cax = cbaxes)
cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist())
cbar.ax.set_ylabel(cbar_label, fontsize = cbar_fontSize)
cbar_tick_labels_pre = [float(t.get_text().replace('−','-')) for t in cbar.ax.get_yticklabels()]
cbar_round = 4
cbar_tick_labels = [round(10**t, cbar_round) for t in cbar_tick_labels_pre]


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



