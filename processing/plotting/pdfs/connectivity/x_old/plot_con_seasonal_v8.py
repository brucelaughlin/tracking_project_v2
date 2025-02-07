# Use the modified pdf files, which have an extra variable for the new box numbers

# v6: combine the progress in v3_1 with the progress in v5

# v3: using .npz files, and make input variable

# V1 copied from ....   Why was I doing transpose before plotting?  That now just feels more confusing.

fig_paramTitle = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall"
#fig_paramTitle = "wc15n model, 300km$^{2}$ coastal boxes, 10km offshore distance as outer wall, physics only, 3D advection, 30-day PLD"

#fig_mainTitle = "Connectivity pdfs (vertical columns integrate to 1).\nSettlement (y-axis) vs. release (x-axis) locations." \
#                "\nGrouped according to season of release."

#fig_mainTitle = "Connectivity pdfs (vertical columns integrate to 1).\nSettlement (y-axis) vs. release (x-axis) locations." \
#                "\nSubplots indicate season of release."

fig_mainTitle = "Connectivity: Log fraction of releases settling.  x-axis = release box, y-axis = settlement box.  Subplots indicate season of release."

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
#fig_fullTitle = fig_mainTitle + "\n" + fig_paramTitle + "\n" + pdf_file_name


base_path = '/home/blaughli/tracking_project/'
pdf_raw_directory = base_path + 'practice/bounding_boxes/final_locations/z_output/z_pre_swap/z_swapped/'

pdf_raw_file = pdf_file_name

d = np.load(pdf_raw_file)

release_counts_per_cell = d['release_counts_per_cell']
pdf_arrays_connectivity = d['pdf_arrays_connectivity']
tick_positions = d['tick_positions']
tick_labels = d['tick_labels']
first_continent_box_dex = d['first_continent_box_num']


# Normalize the histograms along columns (err... rows??) to make connectivity PDFs (for inverse, normalize along rows (err... columns???))

pdf_list = []

pdf_max_val = -999999
pdf_min_val = 999999

for ii in range(len(pdf_arrays_connectivity)):
#for hist in pdf_arrays_connectivity:
    pdf = np.copy(pdf_arrays_connectivity[ii])
    
    pdf = pdf / release_counts_per_cell[ii][:, np.newaxis]

    pdf_list.append(pdf)
    if np.amax(pdf) > pdf_max_val:
        pdf_max_val = np.amax(pdf)
    if np.amin(np.ma.masked_invalid(pdf)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf))

# New approach: set min val to log(0.0001), so 0.0001 is our smallest colorbar tick. indicate that in colorbar tick labels
pdf_min_val = 0.0001

# Also set vmax to log(0.1), to match Patricks's 2011 figures
pdf_max_val = 0.1


# Handling labels and ticks.  Since the box indices are 0-based, and that's how I saved "tick positions", I need to add 1.  

#tick_labels_double = []
tick_labels_double_X = []
tick_labels_double_Y = []

for ii in range(len(tick_labels)):
    tick_labels_double_X.append("{}\n{}".format(tick_positions[ii]+1,tick_labels[ii]))
    if tick_positions[ii]+1 < 10:
        tick_labels_double_Y.append("{}    {}".format(tick_labels[ii],tick_positions[ii]+1))
    else:
        tick_labels_double_Y.append("{}  {}".format(tick_labels[ii],tick_positions[ii]+1))
    if tick_positions[ii]+1 >= first_continent_box_dex:
        tick_positions[ii] = tick_positions[ii] + 1


label_fontsize=6
fig_size = (16,9)

#fig,axs = plt.subplots(2,2, figsize = fig_size, subplot_kw = {'aspect':1})
fig,axs = plt.subplots(2,2, figsize = fig_size)
#plt.setp(axs,xticks=tick_positions,xticklabels=tick_labels_double_X,yticks=tick_positions,yticklabels=tick_labels_double_Y)



num_dummy_lines = 1

ii = 0

# Wait, there are 5 pdfs in the list - the first (index 0) is the overall pdf (non-seasonal).  So, do I have a 1-off error here?

pcs = []

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


    if ii == 1:
        ax = axs[0,0]
        axes_title = "Winter (DJF)"
    elif ii == 2:
        ax = axs[0,1]
        axes_title = "Spring (MAM)"
    elif ii == 3:
        ax = axs[1,0]
        axes_title = "Summer (JJA)"
    else:
        ax = axs[1,1]
        axes_title = "Fall (SON)"

    pcs.append(ax.pcolormesh(X,Y,np.maximum(pdf_separated.T,pdf_min_val),cmap='jet',norm=mpl.colors.LogNorm(vmin=pdf_min_val,vmax=pdf_max_val)))
    ax.plot([0,np.shape(pdf_separated)[1]-1],[0,np.shape(pdf_separated)[0]-1],color="black")
    ax.set_title(axes_title)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels_double_X,fontsize=label_fontsize)
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels_double_Y,fontsize=label_fontsize)






cbar_label = "probability"

cbar = plt.colorbar(pcs[0], ax=axs.ravel(),label=cbar_label,extend="both")

fig.suptitle(fig_fullTitle)

# Create the output directory "figures" if it doesn't exist already
base = os.path.splitext(pdf_raw_file)[0]
figures_directory = base.rsplit('/', 1)[0] + '/figures/'
Path(figures_directory).mkdir(parents=True, exist_ok=True)

fig_file = figures_directory + base.split('/')[-1]  + ".png"

plt.savefig(fig_file)


#plt.show()



