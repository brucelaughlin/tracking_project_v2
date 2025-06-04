# Use the modified pdf files, which have an extra variable for the new box numbers

# V1 copied from ....   Why was I doing transpose before plotting?  That now just feels more confusing.

#---------------------------------------------------------------------
# PDF file - contains labels
#---------------------------------------------------------------------
#pdf_file_name_pre = 'pdf_data_output_releaseLoc_vs_settleTime_test3.p'
pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_test4_physics_only_AKs_1en5.p"

#---------------------------------------------------------------------

#---------------------------------------------------------------------
pdf_file_name = pdf_file_name_pre[0:-2] + "_swapped.p"
#---------------------------------------------------------------------

# Title
#---------------------------------------------------------------------


fig_paramTitle = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"
fig_mainTitle = "Connectivity pdfs (vertical columns integrate to 1).\nSettlement (y-axis) vs. release (x-axis) locations." \
                "\nGrouped according to season of release."
fig_fullTitle = fig_mainTitle + "\n" + fig_paramTitle

cbar_label = "Log base 10 value of probability"
cbar_fontSize = 15
cbar_nBins = 25
#---------------------------------------------------------------------


import pickle
import numpy as np
import matplotlib.pyplot as plt

base_path = '/home/blaughli/tracking_project/'
pdf_raw_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/'

pdf_raw_file = pdf_raw_directory + pdf_file_name
#pdf_raw_file = pdf_raw_directory + 'pdf_data_output_releaseLoc_vs_settleTime_test3.p'




file = open(pdf_raw_file,'rb')
hist_list_exposure_T_source_swapped,hist_list_of_lists_O2_source_swapped,hist_list_connectivity_swapped,hist_list_settleTime_swapped,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels,first_continent_box_dex,oxygen_limit_list = pickle.load(file)
#hist_list_connectivity,hist_list_settleTime,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels = pickle.load(file)
file.close()


#pdf_list = []
#for pdf in pdf_list_connectivity:
#    pdf_list.append(np.log10(np.transpose(pdf)))


# Normalize the histograms along columns (err... rows??) to make connectivity PDFs (for inverse, normalize along rows (err... columns???))

pdf_list = []

pdf_max_val = -999999
pdf_min_val = 999999

for hist in hist_list_connectivity_swapped:
#for hist in hist_list_connectivity:
    pdf = hist
    row_sums = pdf.sum(axis=1)
    pdf = pdf / row_sums[:, np.newaxis]
    pdf = np.log10(pdf)
    #pdf = np.log10(np.transpose(pdf))
    pdf_list.append(pdf)
    if np.amax(pdf) > pdf_max_val:
        pdf_max_val = np.amax(pdf)
    if np.amin(np.ma.masked_invalid(pdf)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf))

print(pdf_max_val)
print(pdf_min_val)

# Determined elsewhere (see/run "check_box_numbers.py")
#first_continent_box_dex = 20
num_dummy_lines = 1


stagger_dex = 0
tick_labels_double_X = []
for ii in range(len(tick_labels)):
    stagger_dex += 1
    if (tick_positions[ii]+1 < first_continent_box_dex) and (stagger_dex % 2 == 0):
    #if (tick_positions[ii]+1 >= 11) and (tick_positions[ii]+1 <= 17) and (stagger_dex % 2 == 0):
        tick_labels_double_X.append("{}\n\n{}".format(tick_positions[ii]+1,tick_labels[ii]))
    else:
        tick_labels_double_X.append("{}\n{}".format(tick_positions[ii]+1,tick_labels[ii]))
    

stagger_dex = 0
tick_labels_double_Y = []
for ii in range(len(tick_labels)):
    stagger_dex += 1
    if (tick_positions[ii]+1 < first_continent_box_dex) and (stagger_dex % 2 == 0):
    #if (tick_positions[ii]+1 >= 11) and (tick_positions[ii]+1 <= 17) and (stagger_dex % 2 == 0):
        tick_labels_double_Y.append("{}       {}".format(tick_labels[ii],tick_positions[ii]+1))
    else:
        tick_labels_double_Y.append("{} {}".format(tick_labels[ii],tick_positions[ii]+1))
    

fig,axs = plt.subplots(2,2)
plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels_double_X,yticks=tick_positions,yticklabels=tick_labels_double_Y)


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
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,0].title.set_text("Winter (DJF)")
        #axs[0,0].set_xlabel("release cell")
        #axs[0,0].set_ylabel("settling cell")
    elif ii == 2:
        mesh2 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,1].title.set_text("Spring (MAM)")
    elif ii == 3:
        mesh3 = axs[1,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,0].title.set_text("Summer (JJA)")
    else:
        mesh4 = axs[1,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,1].title.set_text("Fall (SON)")


#cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist())
#cbar.ax.set_ylabel(cbar_label, fontsize = cbar_fontSize)
#cbar.ax.locator_params(nbins=cbar_nBins)


cbar_label = "Log base 10 of probability"
cbar_fontSize = 15
cbar_nBins = 25

cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist())
cbar.ax.set_ylabel(cbar_label, fontsize = cbar_fontSize)
cbar.ax.locator_params(nbins=cbar_nBins)

cbar_tick_labels_pre = [float(t.get_text().replace('−','-')) for t in cbar.ax.get_yticklabels()]
cbar_round = 4
cbar_tick_labels = [round(10**t, cbar_round) for t in cbar_tick_labels_pre]

cbar_label_2 = "probability"

# define functions that relate the two colorbar scales
# e.g., Celcius to Fahrenheit and vice versa
def logP_to_P(x):
    val = 10**(np.ma.masked_invalid(x))
    return np.ma.masked_invalid(val)
def P_to_logP(x):
    val = np.log10(np.ma.masked_invalid(x))
    return np.ma.masked_invalid(val)

# create a second axes
cbar2 = cbar.ax.secondary_yaxis('left',functions=(logP_to_P,P_to_logP))
cbar2.set_ylabel(cbar_label_2, fontsize = cbar_fontSize, y=0.87)

#cbar_nBins_2 = 15
cbar_nBins_2 = 10

cbar2.locator_params(nbins=cbar_nBins_2)





fig.suptitle(fig_fullTitle)

plt.show()



