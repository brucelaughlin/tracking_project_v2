# Use the modified pdf files, which have an extra variable for the new box numbers

# v5: add 1-1 line

# v4: save figures, don't show.  eventually process entire directory

# v3: using .npz files, and make input variable

# V1 copied from ....   Why was I doing transpose before plotting?  That now just feels more confusing.

fig_paramTitle = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall"
#fig_paramTitle = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"

#fig_mainTitle = "Connectivity pdfs (vertical columns integrate to 1).\nSettlement (y-axis) vs. release (x-axis) locations." \
#                "\nGrouped according to season of release."

#fig_mainTitle = "Connectivity pdfs (vertical columns integrate to 1).\nSettlement (y-axis) vs. release (x-axis) locations." \
#                "\nSubplots indicate season of release."

fig_mainTitle = "Connectivity pdfs (vertical columns integrate to 1).  x-axis = release box, y-axis = settlement box.  Subplots indicate season of release."

#fig_fullTitle = fig_mainTitle + "\n" + fig_paramTitle

#---------------------------------------------------------------------


import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as tkr
import argparse
import os
from pathlib import Path

#---------------------------------------------------------------------
# PDF file - contains labels
#---------------------------------------------------------------------
#pdf_file_name_pre = "pdf_data_output_seasonal_rangeO2_v4_test4_physics_only_AKs_1en5.p"

parser = argparse.ArgumentParser()
parser.add_argument("settlefile")
#parser.add_argument("--settlefile")
args = parser.parse_args()
pdf_file_name = args.settlefile
#pdf_file_name_pre = args.settlefile

#pdf_file_name = pdf_file_name_pre[0:-2] + "_swapped.p"

fig_fullTitle = fig_mainTitle + "\n" + fig_paramTitle + "\n" + os.path.splitext(pdf_file_name.split('/')[-1])[0]
#fig_fullTitle = fig_mainTitle + "\n" + fig_paramTitle + "\n" + pdf_file_name.split('/')[-1]
#fig_fullTitle = fig_mainTitle + "\n" + fig_paramTitle + "\n" + pdf_file_name


base_path = '/home/blaughli/tracking_project/'
pdf_raw_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/z_swapped/'

#pdf_raw_file = pdf_raw_directory + pdf_file_name

pdf_raw_file = pdf_file_name


d = np.load(pdf_raw_file)

pdf_arrays_connectivity = d['pdf_arrays_connectivity']
tick_positions = d['tick_positions']
tick_labels = d['tick_labels']
first_continent_box_dex = d['first_continent_box_num']

#file = open(pdf_raw_file,'rb')
#hist_list_exposure_T_source_swapped,hist_list_of_lists_O2_source_swapped,pdf_arrays_connectivity,hist_list_settleTime_swapped,settlement_boxes_test_array,settlement_times_test_array,counter_array,box_num_mod,tick_positions,tick_labels,first_continent_box_dex,oxygen_limit_list = pickle.load(file)
#file.close()





# Normalize the histograms along columns (err... rows??) to make connectivity PDFs (for inverse, normalize along rows (err... columns???))

pdf_list = []

pdf_max_val = -999999
pdf_min_val = 999999

for hist in pdf_arrays_connectivity:
    pdf = np.copy(hist)
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
    


fig_size = (16,9)

fig,axs = plt.subplots(2,2, figsize = fig_size)
#fig,axs = plt.subplots(2,2)
plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels_double_X,yticks=tick_positions,yticklabels=tick_labels_double_Y)

    


ii = 0

# Wait, there are 5 pdfs in the list - the first (index 0) is the overall pdf (non-seasonal).  So, do I have a 1-off error here?

for pdf_plot in pdf_list[1:]:
#for pdf_plot in pdf_list:

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

#    if ii == 1:
#        axes_obj = axs[0,0]
#        axes_title = "Winter (DJF)"
#    if ii == 2:
#        axes_obj = axs[0,1]
#        axes_title = "Spring (MAM)"
#    if ii == 3:
#        axes_obj = axs[1,0]
#        axes_title = "Summer (JJA)"
#    else:
#        axes_obj = axs[1,1]
#        axes_title = "Fall (SON)"
#
#    mesh1 = axes_obj.pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
#    axes_obj.plot([0,np.shape(pdf_separated)[1]],[0,np.shape(pdf_separated)[0]],color="black")
#    axes_obj.title.set_text(axes_title)
#

    if ii == 1:
        mesh1 = axs[0,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,0].plot([0,np.shape(pdf_separated)[1]],[0,np.shape(pdf_separated)[0]],color="black")
        axs[0,0].title.set_text("Winter (DJF)")
    elif ii == 2:
        mesh2 = axs[0,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[0,1].title.set_text("Spring (MAM)")
    elif ii == 3:
        mesh3 = axs[1,0].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,0].title.set_text("Summer (JJA)")
    else:
        mesh4 = axs[1,1].pcolormesh(X,Y,pdf_separated.T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
        axs[1,1].title.set_text("Fall (SON)")


### hidden plot of overall pdf (not seasonal) for the overall colorbar)
#fig_hide,axs_hide = plt.subplots(1,1, figsize = fig_size)
#mesh_hide = axs_hide.pcolormesh(pdf_list[0].T,cmap='jet',vmin=pdf_min_val,vmax=pdf_max_val)
#fig_hide.set_visible(False)

    



cbar_fontSize = 20

cbar_nBins_2 = 50
#cbar_nBins_2 = 15
#cbar_nBins_2 = 10

last_tick_keep = 10

cbar_label = "probability"

def logP_to_P(x,pos):
    val = round(10**(float(np.ma.masked_invalid(x))),3)
    return val

fmt = mpl.ticker.FuncFormatter(logP_to_P)
#cbar = plt.colorbar(mpl.cm.ScalarMappable(pdf_list[0]), ax=axs.ravel().tolist(), format=fmt)
cbar = plt.colorbar(mesh1, ax=axs.ravel().tolist(), format=fmt)
#cbar = plt.colorbar(mesh_hide, ax=axs.ravel().tolist(), format=fmt)

cbar.ax.set_ylabel(cbar_label, fontsize = cbar_fontSize)
cbar.ax.locator_params(nbins=cbar_nBins_2)
#cbar.ax.yaxis.set_label_position('left')
plt.setp(cbar.ax.get_yticklabels()[0:last_tick_keep], visible=False)

fig.suptitle(fig_fullTitle)

# Create the output directory "figures" if it doesn't exist already

base = os.path.splitext(pdf_raw_file)[0]

figures_directory = base.rsplit('/', 1)[0] + '/figures/'

###figures_directory = pdf_raw_file.split('/')[-2] + '/figures/'

Path(figures_directory).mkdir(parents=True, exist_ok=True)

fig_file = figures_directory + base.split('/')[-1]  + ".png"

#plt.savefig(fig_file)

plt.show()


