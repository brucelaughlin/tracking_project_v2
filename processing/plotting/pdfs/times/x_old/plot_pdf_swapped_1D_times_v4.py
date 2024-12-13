# Copied from "plot_pdf_seasonal_island_swap_v1.py"

# V4: v3 was close, but still not quite

# V3: trying to fix scales between day 1 and rest of days

# V2: try to add a 1-d pdf of day 1 settlement times below 2d plot of days 2-60

#------------------------------------------------------------------
pdf_file_name = "pdf_data_output_releaseLoc_vs_settleTime_test3.p"
#pdf_file_name = "pdf_data_output_seasonal_oneFileTest.p"
#------------------------------------------------------------------


#------------------------------------------------------------------
# Plots are dominated by settlement on day 1 of settlement window - so here's a switch to remove day 1 from the plots
#switch_day_one = True
switch_day_one = False
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

#base_datetime = datetime.datetime(1970,1,1,0,0,0)

file = open(pdf_raw_file,'rb')
pdf_list_connectivity,pdf_list_settleTime,settlement_boxes_test_array,settlement_times_test_array,counter_array = pickle.load(file)  # When the new calc is done, saved a counter_array for checking consistency
file.close()

#pdf_conn_list = []
#for pdf in pdf_list_connectivity:
#    pdf_conn_list.append(np.log10(np.transpose(pdf)))

pdf_max_val = -999999
pdf_min_val = 999999
for pdf in pdf_list_settleTime[1:]:
    pdf_time_full = pdf
    row_sums = pdf_time_full.sum(axis=1)
    pdf_time_full = pdf_time_full / row_sums[:, np.newaxis]
    pdf_time_full = pdf_time_full[:,1:]
    if np.amax(pdf_time_full) > pdf_max_val:
        pdf_max_val = np.amax(pdf_time_full)
    if np.amin(np.ma.masked_invalid(pdf_time_full)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf_time_full))

# Also need min/max for day1 pdfs
pdf_max_val_day1 = -999999
pdf_min_val_day1 = 999999
for pdf in pdf_list_settleTime[1:]:
    pdf_time_full = pdf
    row_sums = pdf_time_full.sum(axis=1)
    pdf_time_full = pdf_time_full / row_sums[:, np.newaxis]
    if np.amax(pdf_time_full) > pdf_max_val_day1:
        pdf_max_val_day1 = np.amax(pdf_time_full)
    if np.amin(np.ma.masked_invalid(pdf_time_full)) < pdf_min_val_day1:
        pdf_min_val_day1 = np.amin(np.ma.masked_invalid(pdf_time_full))




# Determined elsewhere (see/run "check_box_numbers.py")
first_continent_box_dex = 20
num_dummy_lines = 1

tick_positions = [1,5,9,10,13,16,17,19,23,29,32,37,41,47,56,60,67,78]
tick_labels = ['SC','C','SB','SN','SR','A','SC','SM','TJ','PV','PM','PC','PB','CB','PR','PA','CM','CB']

n_boxes_seeded = int(np.shape(pdf_list_settleTime[0])[0]) + num_dummy_lines
n_times = int(np.shape(pdf_list_settleTime[0])[1])
X = np.arange(-0.5, n_boxes_seeded, 1)
if switch_day_one:
    Y = np.arange(-0.5, n_times, 1)
else:
    Y = np.arange(0.5, n_times, 1)

v_scale = 6

y_tick_day1_spacing = .15
y_ticks_day1 = np.arange(pdf_min_val_day1,pdf_max_val_day1, y_tick_day1_spacing)

fig,axs = plt.subplots(4,2, height_ratios = [v_scale,1,v_scale,1], constrained_layout=True)
#plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels,yticks=tick_positions,yticklabels=tick_labels)


for ii in range(1,len(pdf_list_settleTime)):

    pdf_plot_pre1 = pdf_list_settleTime[ii]
    
    row_sums = pdf_plot_pre1.sum(axis=1)
    pdf_plot_pre2 = pdf_plot_pre1 / row_sums[:, np.newaxis]

    # 2D full plot (minus day 1)
    pdf_plot = pdf_plot_pre2[:,1:]
    pdf_separated = np.empty((np.shape(pdf_plot)[0] + num_dummy_lines,np.shape(pdf_plot)[1]))
    pdf_separated[:] = np.nan
    pdf_separated[0:first_continent_box_dex,:] = pdf_plot[0:first_continent_box_dex,:]
    pdf_separated[first_continent_box_dex + num_dummy_lines:,:] = pdf_plot[first_continent_box_dex:,:]

    # 1D day one pdf
    pdf_day1 = pdf_plot_pre2[:,0]
    pdf_day1_separated = np.empty((np.shape(pdf_day1)[0] + num_dummy_lines))
    pdf_day1_separated[:] = np.nan
    pdf_day1_separated[0:first_continent_box_dex] = pdf_day1[0:first_continent_box_dex]
    pdf_day1_separated[first_continent_box_dex + num_dummy_lines:] = pdf_day1[first_continent_box_dex:]

    if ii == 1:
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,0].plot(pdf_day1_separated)
        axs[1,0].set_ylim([pdf_min_val_day1,pdf_max_val_day1])
        axs[1,0].margins(x=0)
        axs[1,0].set_yticks(y_ticks_day1)
        axs[1,0].yaxis.grid(True)
        axs[0,0].title.set_text("Winter (DJF)")
    elif ii == 2:
        mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,1].plot(pdf_day1_separated)
        axs[1,1].set_ylim([pdf_min_val_day1,pdf_max_val_day1])
        axs[1,1].margins(x=0)
        axs[1,1].set_yticks(y_ticks_day1)
        axs[1,1].yaxis.grid(True)
        axs[0,1].title.set_text("Spring (MAM)")
    elif ii == 3:
        mesh1 = axs[2,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[3,0].plot(pdf_day1_separated)
        axs[3,0].set_ylim([pdf_min_val_day1,pdf_max_val_day1])
        axs[3,0].margins(x=0)
        axs[3,0].set_yticks(y_ticks_day1)
        axs[3,0].yaxis.grid(True)
        axs[2,0].title.set_text("Summer (JJA)")
    else:
        mesh1 = axs[2,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[3,1].plot(pdf_day1_separated)
        axs[3,1].set_ylim([pdf_min_val_day1,pdf_max_val_day1])
        axs[3,1].margins(x=0)
        axs[3,1].set_yticks(y_ticks_day1)
        axs[3,1].yaxis.grid(True)
        axs[2,1].title.set_text("Fall (SON)")

    plt.colorbar(mesh1)
if switch_day_one:
    fig.suptitle("PDFs of Settlement Time (y-axis) vs Release Location (x-axis).\n Grouped according to season of release")
else:
    fig.suptitle("PDFs of Settlement Time (y-axis) vs Release Location (x-axis).\n Grouped according to season of release.\n(Day 1 bin removed)")

plt.show()



