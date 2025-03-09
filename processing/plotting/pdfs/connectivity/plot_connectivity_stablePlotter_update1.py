# Copied from plot_con_seasonal_pDrakeCompare_v3.py, which was stable and reproduced Patrick's plots.  
# Many updates I did since then have had issues.

annualOnlySwitch = True

script_version = "stablePlotter_update1"

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
parser.add_argument("pdfRawFile")
#parser.add_argument("--settlefile")
args = parser.parse_args()
pdf_raw_file = args.pdfRawFile

# Create the output directory "figures" if it doesn't exist already
base = os.path.splitext(pdf_raw_file)[0]
figures_directory = base.rsplit('/', 1)[0] + '/figures/'
Path(figures_directory).mkdir(parents=True, exist_ok=True)
file_name_stem = base.split('/')[-1]
fig_file = figures_directory + file_name_stem  + '_' + script_version  + ".png"

wordcount_line_limit = 75
num_lines_filename = int(np.ceil(len(file_name_stem)/wordcount_line_limit))
#num_lines_filename = len(file_name_stem) % wordcount_line_limit
figure_title = ""

for ii in range(num_lines_filename):
    figure_title += f"{file_name_stem[ii*wordcount_line_limit:(ii+1)*wordcount_line_limit]}"
    if ii < num_lines_filename - 1:
        figure_title += " ...\n"

#print(figure_title)

d = np.load(pdf_raw_file)

release_counts_per_cell = d['release_counts_per_cell']
connectivity_histogram_array = d['pdf_arrays_connectivity']

pld = d['pld_days']

print('minimum number of releases in a cell: {}'.format(np.min(release_counts_per_cell)))


num_floats = np.sum(release_counts_per_cell)
num_settlers = np.sum(connectivity_histogram_array)
settle_strength = num_settlers/num_floats

release_months_list = ["Jan-Dec", "Dec-Feb", "March-May","June-Aug","Sep-Nov"]

connectivity_pdf_list = []

pdf_max_val = -999999
pdf_min_val = 999999

for ii in range(len(connectivity_histogram_array)):
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

fig,ax = plt.subplots(1,1, figsize = fig_size)
#fig,axs = plt.subplots(2,2, figsize = fig_size)



num_dummy_lines = 1

# Wait, there are 5 pdfs in the list - the first (index 0) is the overall pdf (non-seasonal).  So, do I have a 1-off error here?

pcs = []


# New approach: set min val to log(0.0001), so 0.0001 is our smallest colorbar tick. indicate that in colorbar tick labels
pdf_min_val = 0
#pdf_min_val = 0.0001



boundary_index = 27
island_index = 489

print(f"Number of pdfs: {len(connectivity_pdf_list)}")

#for pdf_plot in connectivity_pdf_list[1:]:
#for pdf_plot in connectivity_pdf_list:
for pdf_index in range(len(connectivity_pdf_list)):
    
    pdf_plot = connectivity_pdf_list[pdf_index]



    subplot_title = ""
    subplot_title += f"PLD: {pld[0]}-{pld[1]} days"
    subplot_title += f"\nReleased {release_months_list[pdf_index]}, 0-20m"
    subplot_title += f"\nn floats: {num_floats}, n settlers: {num_settlers}, settle strength: {settle_strength:.3f}"

    n_boxes_seeded = int(np.shape(pdf_plot)[1])
    n_boxes_settled = int(np.shape(pdf_plot)[0])
    X = np.arange(-0.5, n_boxes_settled, 1)
    Y = np.arange(-0.5, n_boxes_seeded, 1)

    
    pdf_plot[0:boundary_index,:] = 0
    pdf_plot[:,0:boundary_index] = 0
    pdf_plot[island_index - boundary_index:island_index,:] = 0
    pdf_plot[:,island_index - boundary_index:island_index] = 0

    np.savetxt('connectivity.csv', pdf_plot, delimiter = ',')

    pcs.append(ax.pcolormesh(X,Y,pdf_plot.T,cmap='jet'))
    #pcs.append(ax.pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=mpl.colors.Normalize(vmin=pdf_min_val,vmax=pdf_max_val)))
    
    ax.plot([0,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)

    ax.set_title(subplot_title)
    #ax.title.set(subplot_title)

    if annualOnlySwitch:
        break




cbar_label = "probability"

cbar = plt.colorbar(pcs[0],label=cbar_label,extend="both")
#cbar = plt.colorbar(pcs[0], ax=axs.ravel(),label=cbar_label,extend="both")

ax.set_aspect('equal')

#fig.suptitle(file_name_stem)
fig.suptitle(figure_title)


plt.savefig(fig_file)


#plt.show()



