# Copied from "plot_pdf_seasonal_island_swap_v1.py"

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

for pdf in pdf_list_settleTime:
    if switch_day_one:
        pdf_time_full = pdf
        #pdf_time_full = pdf_list_settleTime[0]
    else:
        pdf_time_full = pdf[:,1:]
        #pdf_time_full = pdf_list_settleTime[0][:,1:]

    row_sums = pdf_time_full.sum(axis=1)
    pdf_time_full = pdf_time_full / row_sums[:, np.newaxis]

    if np.amax(pdf_time_full) > pdf_max_val:
        pdf_max_val = np.amax(pdf_time_full)
    if np.amin(np.ma.masked_invalid(pdf_time_full)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf_time_full))

#print(pdf_max_val)
#print(pdf_min_val)

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

fig,axs = plt.subplots(2,2)
#plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels,yticks=tick_positions,yticklabels=tick_labels)


for ii in range(1,len(pdf_list_settleTime)):
#for pdf_plot in pdf_list:

    if switch_day_one:
        pdf_plot_pre = pdf_list_settleTime[ii]
    else:
        pdf_plot_pre = pdf_list_settleTime[ii][:,1:]

    row_sums = pdf_plot_pre.sum(axis=1)
    pdf_plot = pdf_plot_pre / row_sums[:, np.newaxis]

    pdf_separated = np.empty((np.shape(pdf_plot)[0] + num_dummy_lines,np.shape(pdf_plot)[1]))
    #pdf_separated = np.empty((np.shape(pdf_plot)[0] + num_dummy_lines,np.shape(pdf_plot)[1] + num_dummy_lines))
    pdf_separated[:] = np.nan

    pdf_separated[0:first_continent_box_dex,:] = pdf_plot[0:first_continent_box_dex,:]
    pdf_separated[first_continent_box_dex + num_dummy_lines:,:] = pdf_plot[first_continent_box_dex:,:]

    if ii == 1:
        #mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet')
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,0].title.set_text("Winter (DJF)")
    elif ii == 2:
        #mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet')
        mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,1].title.set_text("Spring (MAM)")
    elif ii == 3:
        #mesh1 = axs[1,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet')
        mesh1 = axs[1,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,0].title.set_text("Summer (JJA)")
    else:
        #mesh1 = axs[1,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet')
        mesh1 = axs[1,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,1].title.set_text("Fall (SON)")

    plt.colorbar(mesh1)
if switch_day_one:
    fig.suptitle("PDFs of Settlement Time (y-axis) vs Release Location (x-axis).\n Grouped according to season of release")
else:
    fig.suptitle("PDFs of Settlement Time (y-axis) vs Release Location (x-axis).\n Grouped according to season of release.\n(Day 1 bin removed)")

plt.show()



