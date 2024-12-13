# Copied from ""

#------------------------------------------------------------------
pdf_file_name_pre = "pdf_data_output_seasonal_v3_oneFileTest.p"

#------------------------------------------------------------------

#------------------------------------------------------------------
#pdf_file_name = pdf_file_name_pre[0:-2] + "_swapped.p"
pdf_file_name = pdf_file_name_pre
#------------------------------------------------------------------


#------------------------------------------------------------------
fig_paramTitle = "wc_15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"
fig_mainTitle = "PDFs of Settlement Time (y-axis) vs Release Location (x-axis).\nGrouped according to season of release."

fig_fullTitle = fig_mainTitle + "\n" + fig_paramTitle

yLabel = "Days since end of PLD"



cbar_label = "Log base 10 value of probability"
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
pdf_list_exposure_T_source,pdf_list_exposure_oxygen_source,pdf_list_connectivity,pdf_list_settleTime,settlement_boxes_test_array,settlement_times_test_array,counter_array = pickle.load(file)
#pdf_list_exposure_T_source,pdf_list_exposure_oxygen_source,pdf_list_connectivity,pdf_list_settleTime,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels = pickle.load(file)
file.close()

stats_vectors_hist = np.zeros((len(pdf_list_exposure_T_source)-1,np.shape(pdf_list_exposure_T_source[0])[1]))
stats_vectors = np.zeros((len(pdf_list_exposure_T_source)-1,np.shape(pdf_list_exposure_T_source[0])[1]))

pdf_list_full = []

pdf_max_val = -999999
pdf_min_val = 999999
ii = 0
for pdf in pdf_list_exposure_T_source[1:]:
    pdf_time_full = pdf

    stats_vectors_hist[ii,:] = pdf_time_full.sum(axis=0)
    ii += 1

    row_sums = pdf_time_full.sum(axis=1)
    pdf_time_full = pdf_time_full / row_sums[:, np.newaxis]
    pdf_time_full = np.log10(pdf_time_full)
    pdf_list_full.append(pdf_time_full)
    if np.amax(pdf_time_full) > pdf_max_val:
        pdf_max_val = np.amax(pdf_time_full)
    if np.amin(np.ma.masked_invalid(pdf_time_full)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf_time_full))

# Calc the stats
for ii in range(np.shape(stats_vectors_hist)[0]):
    stats_vectors[ii,:] = stats_vectors_hist[ii,:] / np.sum(stats_vectors_hist[ii,:])
    stats_vectors[ii,:] = np.cumsum(stats_vectors[ii,:])


#stats_10 = list(range(int(np.shape(stats_vectors)[0])))
#stats_90 = list(range(int(np.shape(stats_vectors)[0])))

stats_10 = []
stats_50 = []
stats_60 = []
stats_90 = []

compare_val = .10
for ii in range(np.shape(stats_vectors)[0]):
    for jj in range(np.shape(stats_vectors)[1]):
        if stats_vectors[ii,jj] > compare_val:
            stats_10.append(jj)
            break

compare_val = .50
for ii in range(np.shape(stats_vectors)[0]):
    for jj in range(np.shape(stats_vectors)[1]):
        if stats_vectors[ii,jj] > compare_val:
            stats_50.append(jj)
            break

compare_val = .60
for ii in range(np.shape(stats_vectors)[0]):
    for jj in range(np.shape(stats_vectors)[1]):
        if stats_vectors[ii,jj] > compare_val:
            stats_60.append(jj)
            break

compare_val = .90
for ii in range(np.shape(stats_vectors)[0]):
    for jj in range(np.shape(stats_vectors)[1]):
        if stats_vectors[ii,jj] > compare_val:
            stats_90.append(jj)
            break


# Determined elsewhere (see/run "check_box_numbers.py")
first_continent_box_dex = 20
num_dummy_lines = 1

n_boxes_seeded = int(np.shape(pdf_list_exposure_T_source[0])[0]) + num_dummy_lines
n_times = int(np.shape(pdf_list_exposure_T_source[0])[1])
X = np.arange(-0.5, n_boxes_seeded, 1)
Y = np.arange(-0.5, n_times, 1)

v_scale = 6



stagger_dex = 0
tick_labels_double_X = []
for ii in range(len(tick_labels)):
    stagger_dex += 1
    if (tick_positions[ii]+1 >= 11) and (tick_positions[ii]+1 <= 17) and (stagger_dex % 2 == 0):
        tick_labels_double_X.append("{}\n\n{}".format(tick_positions[ii]+1,tick_labels[ii]))
    else:
        tick_labels_double_X.append("{}\n{}".format(tick_positions[ii]+1,tick_labels[ii]))


fig,axs = plt.subplots(2,2)

plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels_double_X)

label_fontSize = 10
#label_fontSize = 8

for ii in range(len(pdf_list_full)):
#for ii in range(1,len(pdf_list_exposure_T_source)):

    #pdf_plot_pre1 = pdf_list_exposure_T_source[ii]
    
    #row_sums = pdf_plot_pre1.sum(axis=1)
    #pdf_plot_pre2 = pdf_plot_pre1 / row_sums[:, np.newaxis]

    # 2D full plot (minus day 1)
    pdf_plot = pdf_list_full[ii]
    #pdf_plot = pdf_plot_pre2[:,1:]
    pdf_separated = np.empty((np.shape(pdf_plot)[0] + num_dummy_lines,np.shape(pdf_plot)[1]))
    pdf_separated[:] = np.nan
    pdf_separated[0:first_continent_box_dex,:] = pdf_plot[0:first_continent_box_dex,:]
    pdf_separated[first_continent_box_dex + num_dummy_lines:,:] = pdf_plot[first_continent_box_dex:,:]

    if ii == 1:
        #mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        pname = "Winter (DJF)"
        ptitle = "50% by day {}; 90% by day {};    {}".format(stats_50[ii]+1,stats_90[ii]+1,pname)
        axs[0,0].title.set_text(ptitle)
        #axs[0,0].title.set_text("Winter (DJF)")
        axs[0,0].set_ylabel(yLabel)
        axs[0,0].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        axs[0,0].set_xticks(tick_positions)
        #axs[0,0].set_xticklabels(tick_labels, fontsize=label_fontSize)
        axs[0,0].set_xticklabels(tick_labels_double_X, fontsize=label_fontSize)
        #------------------------------------------------------------------
    elif ii == 2:
        #mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        mesh1 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        pname = "Spring (MAM)"
        ptitle = "50% by day {}; 90% by day {};    {}".format(stats_50[ii]+1,stats_90[ii]+1,pname)
        axs[0,1].title.set_text(ptitle)
        #axs[0,1].title.set_text("Spring (MAM)")
        axs[0,1].set_ylabel(yLabel)
        axs[0,1].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        axs[0,1].set_xticks(tick_positions)
        #axs[0,1].set_xticklabels(tick_labels, fontsize=label_fontSize)
        axs[0,1].set_xticklabels(tick_labels_double_X, fontsize=label_fontSize)
        #------------------------------------------------------------------
    elif ii == 3:
        #mesh1 = axs[1,0].pcolormesh(X,Y,pdf_separated,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        mesh1 = axs[1,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        pname = "Summer (JJA)"
        ptitle = "50% by day {}; 90% by day {};    {}".format(stats_50[ii]+1,stats_90[ii]+1,pname)
        axs[1,0].title.set_text(ptitle)
        #axs[1,0].title.set_text("Summer (JJA)")
        axs[1,0].set_ylabel(yLabel)
        axs[1,0].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        axs[1,0].set_xticks(tick_positions)
        #axs[1,0].set_xticklabels(tick_labels, fontsize=label_fontSize)
        axs[1,0].set_xticklabels(tick_labels_double_X, fontsize=label_fontSize)
        #------------------------------------------------------------------
    else:
        #mesh1 = axs[1,1].pcolormesh(X,Y,pdf_separated,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        mesh1 = axs[1,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        pname = "Fall (SON)"
        ptitle = "50% by day {}; 90% by day {};    {}".format(stats_50[ii]+1,stats_90[ii]+1,pname)
        axs[1,1].title.set_text(ptitle)
        #axs[1,1].title.set_text("Fall (SON)")
        axs[1,1].set_ylabel(yLabel)
        axs[1,1].yaxis.label.set(fontsize=15)
        #------------------------------------------------------------------
        axs[1,1].set_xticks(tick_positions)
        #axs[1,1].set_xticklabels(tick_labels, fontsize=label_fontSize)
        axs[1,1].set_xticklabels(tick_labels_double_X, fontsize=label_fontSize)
        #------------------------------------------------------------------

#stats_10 = list(range(int(np.shape(stats_vectors)[0])))
#stats_90 = list(range(int(np.shape(stats_vectors)[0])))

    #plt.colorbar(mesh1)

#fig.colorbar(mesh1, ax=axs.ravel().tolist())

cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist())
cbar.ax.set_ylabel(cbar_label, fontsize = cbar_fontSize)
cbar.ax.locator_params(nbins=cbar_nBins)


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



