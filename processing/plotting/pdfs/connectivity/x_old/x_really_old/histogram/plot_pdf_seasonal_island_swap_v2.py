# Use the modified pdf files, which have an extra variable for the new box numbers

# V2: add 2nd layer of axis lables, for island names

import pickle
import numpy as np
import matplotlib.pyplot as plt

base_path = '/home/blaughli/tracking_project/'
pdf_raw_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/'

pdf_raw_file = pdf_raw_directory + 'pdf_data_output_releaseLoc_vs_settleTime_test3.p'
#pdf_raw_file = pdf_raw_directory + 'pdf_data_output_seasonal_test3_modified.p'
#pdf_raw_file = pdf_raw_directory + 'pdf_data_output_seasonal_test3.p'
#pdf_raw_file = pdf_raw_directory + 'pdf_data_output_originalOrder.p'
#pdf_raw_file = pdf_raw_directory + 'pdf_data_output.p'
#pdf_raw_file = pdf_raw_directory + 'pdf_data_output_seasonal.p'
#pdf_raw_file = pdf_raw_directory + 'test_pdf_data.p'




file = open(pdf_raw_file,'rb')
pdf_list_connectivity,pdf_list_settleTime,settlement_boxes_test_array,settlement_times_test_array,counter_array = pickle.load(file)  # When the new calc is done, saved a counter_array for checking consistency
#pdf_raw,pdf_raw_djf,pdf_raw_mam,pdf_raw_jja,pdf_raw_son,counter_array,box_num_islands_mod = pickle.load(file)  # When the new calc is done, saved a counter_array for checking consistency
#pdf_raw,pdf_raw_djf,pdf_raw_mam,pdf_raw_jja,pdf_raw_son,counter_array = pickle.load(file)  # When the new calc is done, saved a counter_array for checking consistency
file.close()





pdf_list = []
for pdf in pdf_list_connectivity:
    pdf_list.append(np.log10(np.transpose(pdf)))


pdf_max_val = -999999
pdf_min_val = 999999
for pdf in pdf_list:
    if np.amax(pdf) > pdf_max_val:
        pdf_max_val = np.amax(pdf)
    if np.amin(np.ma.masked_invalid(pdf)) < pdf_min_val:
    #if (np.amin(pdf) < pdf_min_val and np.amin(pdf) >= 0):
    #if np.amin(pdf) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf))
        #pdf_min_val = np.amin(pdf)

print(pdf_max_val)
print(pdf_min_val)

# Determined elsewhere (see/run "check_box_numbers.py")
first_continent_box_dex = 20
num_dummy_lines = 1

#tick_positions = [0,4,8,10,12,23,29,32,37,41,47,56,60,67,78]
#tick_labels = ['SC','C','SB','SN','SM','TJ','PV','PM','PC','PB','CB','PR','PA','CM','CB']
tick_positions = [1,5,9,10,13,16,17,19,23,29,32,37,41,47,56,60,67,78]
tick_labels = ['SC','C','SB','SN','SR','A','SC','SM','TJ','PV','PM','PC','PB','CB','PR','PA','CM','CB']

tick_labels_double_X = []
for ii in range(len(tick_labels)):
    tick_labels_double_X.append("{}\n{}".format(tick_positions[ii],tick_labels[ii]))

tick_labels_double_Y = []
for ii in range(len(tick_labels)):
    tick_labels_double_Y.append("{} {}".format(tick_labels[ii],tick_positions[ii]))

fig,axs = plt.subplots(2,2)
plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels_double_X,yticks=tick_positions,yticklabels=tick_labels_double_Y)
#plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels_double,yticks=tick_positions,yticklabels=tick_labels)
#plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels,yticks=tick_positions,yticklabels=tick_labels)


ii = 0
for pdf_plot in pdf_list:

    ii += 1

    pdf_separated = np.empty((np.shape(pdf_plot)[0] + num_dummy_lines,np.shape(pdf_plot)[1] + num_dummy_lines))
    pdf_separated[:] = np.nan

    pdf_separated[0:first_continent_box_dex,0:first_continent_box_dex] = pdf_plot[0:first_continent_box_dex,0:first_continent_box_dex]
    pdf_separated[first_continent_box_dex + num_dummy_lines:,0:first_continent_box_dex] = pdf_plot[first_continent_box_dex:,0:first_continent_box_dex]
    pdf_separated[0:first_continent_box_dex,first_continent_box_dex + num_dummy_lines:] = pdf_plot[0:first_continent_box_dex,first_continent_box_dex:]
    pdf_separated[first_continent_box_dex + num_dummy_lines:,first_continent_box_dex + num_dummy_lines:] = pdf_plot[first_continent_box_dex:,first_continent_box_dex:]

    n_boxes_seeded = int(np.shape(pdf_separated)[1])
    n_boxes_settled = int(np.shape(pdf_separated)[0])
    X = np.arange(-0.5, n_boxes_settled, 1)
    Y = np.arange(-0.5, n_boxes_seeded, 1)

    if ii == 1:
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,0].title.set_text("Winter (DJF)")
    elif ii == 2:
        mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,1].title.set_text("Spring (MAM)")
    elif ii == 3:
        mesh1 = axs[1,0].pcolormesh(X,Y,pdf_separated,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,0].title.set_text("Summer (JJA)")
    else:
        mesh1 = axs[1,1].pcolormesh(X,Y,pdf_separated,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,1].title.set_text("Fall (SON)")

    plt.colorbar(mesh1)


fig.suptitle("PDFs of Settlement (y-axis) vs Release (x-axis) locations.\n Grouped according to season of release")

plt.show()



