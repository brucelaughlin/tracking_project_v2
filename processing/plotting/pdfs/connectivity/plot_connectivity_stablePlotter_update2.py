# Copied from plot_con_seasonal_pDrakeCompare_v3.py, which was stable and reproduced Patrick's plots.  
# Many updates I did since then have had issues.

# update2: add input argument to specify whether using patrick cells or not.  if not, don't do the trimming done in versions
# assuming use of patrick cells

annualOnlySwitch = True

script_version = "stablePlotter_update2"

# Also set vmax to log(0.1), to match Patricks's 2011 figures
#pdf_max_val = 0.02
#pdf_max_val = 0.01
#pdf_max_val = 0.1
#pdf_max_val = 0.2
#---------------------------------------------------------------------


import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as tkr
import argparse
import os
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("connectivitydatafile")
parser.add_argument("pdrakeswitch", nargs='?', type=str)
args = parser.parse_args()

connectivity_data_file = args.connectivitydatafile
pDrake_switch = args.pdrakeswitch

output_csv_dir_stem = "y_csv_files"
output_figures_dir_stem = "z_figures"

connectivity_data_dir = os.path.dirname(os.path.realpath(connectivity_data_file))
connectivity_file_stem = Path(connectivity_data_file).stem


# Create the output directories if they don't exist already
#base = os.path.splitext(connectivity_data_file)[0]

output_figures_dir = os.path.join(connectivity_data_dir,output_figures_dir_stem)
#output_figures_dir = base.rsplit('/', 1)[0] + "/" + output_figures_dir_stem
Path(output_figures_dir).mkdir(parents=True, exist_ok=True)

output_csv_dir = os.path.join(connectivity_data_dir,output_csv_dir_stem)
#output_csv_dir = base.rsplit('/', 1)[0] + "/" + output_csv_dir_stem
Path(output_csv_dir).mkdir(parents=True, exist_ok=True)


figure_file_leaf = connectivity_file_stem + '.png'
figure_file = os.path.join(output_figures_dir,figure_file_leaf)

csv_file_leaf = connectivity_file_stem + '.png'
csv_file = os.path.join(output_csv_dir,csv_file_leaf)


wordcount_line_limit = 75
num_lines_filename = int(np.ceil(len(connectivity_file_stem)/wordcount_line_limit))
#num_lines_filename = int(np.ceil(len(file_name_stem)/wordcount_line_limit))
#num_lines_filename = len(file_name_stem) % wordcount_line_limit
figure_title = ""

for ii in range(num_lines_filename):
    figure_title += f"{connectivity_file_stem[ii*wordcount_line_limit:(ii+1)*wordcount_line_limit]}"
    #figure_title += f"{file_name_stem[ii*wordcount_line_limit:(ii+1)*wordcount_line_limit]}"
    if ii < num_lines_filename - 1:
        figure_title += " ...\n"

#print(figure_title)

d = np.load(connectivity_data_file)

release_counts_per_cell = d['release_counts_per_cell']
connectivity_histogram_array = d['pdf_arrays_connectivity']

num_pdfs = len(connectivity_histogram_array)

pld = d['pld_days']

print('minimum number of releases in a cell: {}'.format(np.min(release_counts_per_cell)))


if annualOnlySwitch:
    pdfs_to_plot_indices = [0]
else:
    pdfs_to_plot_indices = [1,2,3,4]


num_floats = int(np.sum(release_counts_per_cell))
num_settlers = int(np.sum(connectivity_histogram_array))
settle_strength = num_settlers/num_floats

release_months_list = ["Jan-Dec", "Dec-Feb", "March-May","June-Aug","Sep-Nov"]

connectivity_pdf_list = []

pdf_max_val = -999999
pdf_min_val = 999999

for ii in pdfs_to_plot_indices:
#for ii in range(len(connectivity_histogram_array)):
    pdf = np.copy(connectivity_histogram_array[ii])
   
    #print(np.count_nonzero(np.isnan(pdf)))
    
    #for jj in range(np.shape(pdf)[0]):
        #if pdf[jj,jj] > 1:
        #    print(f'{jj}: {pdf[jj,jj]}')


#    import pdb; pdb.set_trace()

    # recall paul's lesson; np.newaxis really means None, and this is about the pointwise expansion business... very convenient
    pdf = pdf / release_counts_per_cell[ii][:, np.newaxis]

    #print(np.count_nonzero(np.isnan(pdf)))
    
    # trim the pdfs, to match Patrick's plots..
    #pdf = pdf[cutoff_dex_tj:,cutoff_dex_tj:]
    #pdf = pdf[cutoff_dex_pv:,cutoff_dex_tj:]

    connectivity_pdf_list.append(pdf)
    if np.amax(pdf) > pdf_max_val:
        pdf_max_val = np.amax(pdf)
    if np.amin(np.ma.masked_invalid(pdf)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf))



label_fontsize=6
fig_size = (16,9)

if annualOnlySwitch:
    fig,ax = plt.subplots(1,1, figsize = fig_size)
else:
    fig,ax = plt.subplots(2,2, figsize = fig_size)



num_dummy_lines = 1

# Wait, there are 5 pdfs in the list - the first (index 0) is the overall pdf (non-seasonal).  So, do I have a 1-off error here?

pcs = []


# New approach: set min val to log(0.0001), so 0.0001 is our smallest colorbar tick. indicate that in colorbar tick labels
pdf_min_val = 0
#pdf_min_val = 0.0001

pdf_max_val = 0.01

boundary_index = 27
island_index = 489

print(f"Number of pdfs: {len(connectivity_pdf_list)}")

#for pdf_plot in connectivity_pdf_list[1:]:
#for pdf_plot in connectivity_pdf_list:
#for pdf_index in range(len(connectivity_pdf_list)):
for pdf_index in pdfs_to_plot_indices:
   
    if pdf_index > 0:
        axis_index = pdf_index - 1
    else:
        axis_index = 0

    pdf_plot = connectivity_pdf_list[pdf_index]



    subplot_title = ""
    subplot_title += f"PLD: {pld[0]}-{pld[1]} days"
    subplot_title += f"\nReleased {release_months_list[pdf_index]}, 0-20m"
    subplot_title += f"\nn floats: {num_floats}, n settlers: {num_settlers}, settle strength: {settle_strength:.3f}"

    n_boxes_seeded = int(np.shape(pdf_plot)[1])
    n_boxes_settled = int(np.shape(pdf_plot)[0])
    X = np.arange(-0.5, n_boxes_settled, 1)
    Y = np.arange(-0.5, n_boxes_seeded, 1)

    #if pDrake_switch is None:
    #    pdf_plot[0:boundary_index,:] = 0
    #    pdf_plot[:,0:boundary_index] = 0
    #    pdf_plot[island_index - boundary_index:island_index,:] = 0
    #    pdf_plot[:,island_index - boundary_index:island_index] = 0

    if pdf_index == 0:
        conn_csv_file_name = output_csv_dir_stem + "/" + connectivity_file_stem + f"_annual.csv"
        #conn_csv_file_name = output_csv_dir_stem + "/" + file_name_stem + f"_connectivityData_annual.csv"
    else:
        conn_csv_file_name = output_csv_dir_stem + "/" + connetivity_file_stem + f"_season_{pdf_index + 1}/4.csv"
        #conn_csv_file_name = output_csv_dir_stem + "/" + file_name_stem + f"_connectivityData_season_{pdf_index + 1}/4.csv"

    ###np.savetxt('connectivity.csv', pdf_plot, delimiter = ',')
    #np.savetxt(conn_csv_file_name, pdf_plot, delimiter = ',')

    #pcs.append(ax.pcolormesh(X,Y,pdf_plot.T,cmap='jet'))
    pcs.append(ax.pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=mpl.colors.Normalize(vmin=pdf_min_val,vmax=pdf_max_val)))
    
    ax.plot([axis_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)
    #ax.plot([0,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)

    ax.set_title(subplot_title)
    #ax.title.set(subplot_title)

    #if annualOnlySwitch:
    #    break




cbar_label = "probability"

cbar = plt.colorbar(pcs[0],label=cbar_label,extend="both")
#cbar = plt.colorbar(pcs[0], ax=axs.ravel(),label=cbar_label,extend="both")

ax.set_aspect('equal')

#fig.suptitle(figure_title)

plt.figtext(0.5, 0.02, figure_title, ha="center")

plt.savefig(figure_file)


#plt.show()



