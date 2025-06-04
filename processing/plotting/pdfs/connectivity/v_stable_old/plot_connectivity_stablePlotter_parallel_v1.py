# Copied from "stablePlotter_update2", modified to handle multiple input files from the parallelized connectivity histogram calculation



# ----------------------------------------------
# ----------------------------------------------
# Only for testing
# ----------------------------------------------
# ----------------------------------------------

annualOnlySwitch = True

# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------

script_version = "stablePlotter_parallel_v1"

# Also set vmax to log(0.1), to match Patricks's 2011 figures
#pdf_max_val = 0.02
#pdf_max_val = 0.01
#pdf_max_val = 0.1
#pdf_max_val = 0.2
#---------------------------------------------------------------------

import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("trackingdir")
#parser.add_argument("connectivitydatadir")
parser.add_argument("pdrakeswitch", nargs='?', type=str)
args = parser.parse_args()

tracking_dir = args.trackingdir
#connectivity_data_dir = args.connectivitydatadir
pDrake_switch = args.pdrakeswitch

input_connectivity_dir_stem = "v_connHist_files"
output_csv_dir_stem = "y_csv_files"
output_figures_dir_stem = "z_figures"

connectivity_data_dir = os.path.join(tracking_dir,input_connectivity_dir_stem)


output_figures_dir = os.path.join(connectivity_data_dir,output_figures_dir_stem)
Path(output_figures_dir).mkdir(parents=True, exist_ok=True)

output_csv_dir = os.path.join(connectivity_data_dir,output_csv_dir_stem)
Path(output_csv_dir).mkdir(parents=True, exist_ok=True)

figure_title_base = Path(tracking_dir).stem

connectivity_pdf_list = []
pdf_max_val = -999999
pdf_min_val = 999999
release_months_list = ["Jan-Dec", "Dec-Feb", "March-May","June-Aug","Sep-Nov"]
    
if annualOnlySwitch:
    pdfs_to_plot_indices = [0]
else:
    pdfs_to_plot_indices = [1,2,3,4]



num_floats = 0
num_settlers = 0


# Copying some SO code to loop over files
directory = os.fsencode(connectivity_data_dir)

file_counter = 0

for connectivity_data_file_byte in os.listdir(directory):

    connectivity_data_file = os.fsdecode(connectivity_data_file_byte)
    connectivity_data_dir = os.fsdecode(directory)

    connectivity_data_file = os.path.join(connectivity_data_dir,connectivity_data_file)

    if connectivity_data_file.endswith(".npz"):

#        print('hi4')

        file_counter += 1


#        print(f"FILE: {connectivity_data_file}")
#        print(f"file counter: {file_counter}")

        if file_counter == 1:

            connectivity_file_stem = Path(connectivity_data_file).stem

            # The figure titles also need to change when I start plotting seasons again.  So, this is bugged
            # for seasons

            figure_file_leaf = figure_title_base + '.png'
            #figure_file_leaf = connectivity_file_stem + '.png'
            figure_file = os.path.join(output_figures_dir,figure_file_leaf)

            csv_file_leaf = connectivity_file_stem + '.png'
            csv_file = os.path.join(output_csv_dir,csv_file_leaf)

            wordcount_line_limit = 75
            num_lines_filename = int(np.ceil(len(connectivity_file_stem)/wordcount_line_limit))
            figure_title = ""
            
            # This code will not really work for seasonal plots... been debugging with the annual pdf for so long, I've
            # written bugs into this script.

            for ii in range(num_lines_filename):
                figure_title += f"{figure_title_base[ii*wordcount_line_limit:(ii+1)*wordcount_line_limit]}"
                #figure_title += f"{connectivity_file_stem[ii*wordcount_line_limit:(ii+1)*wordcount_line_limit]}"
                if ii < num_lines_filename - 1:
                    figure_title += " ...\n"

            # Setup stuff, now that we're combining multiple histograms 
            d = np.load(connectivity_data_file)
            pld = d['pld_days']
            release_counts_per_cell = d['release_counts_per_cell']
            connectivity_histogram_array = d['pdf_arrays_connectivity']
            num_pdfs = len(connectivity_histogram_array)
            
            release_counts_per_cell_list = []
 #           print('hi2')
            for ii in pdfs_to_plot_indices:
 #               print('hi3')
                pdf = np.copy(connectivity_histogram_array[ii])
                connectivity_pdf_list.append(np.zeros_like(pdf))
                release_counts_per_cell_list.append(np.zeros_like(release_counts_per_cell))
               

        
        # Re-loading on the first loop iteration seems sloppy, but I don't see the harm 
        d = np.load(connectivity_data_file)

        release_counts_per_cell = d['release_counts_per_cell']
        connectivity_histogram_array = d['pdf_arrays_connectivity']

        #print('minimum number of releases in a cell: {}'.format(np.min(release_counts_per_cell)))

        num_floats += int(np.sum(release_counts_per_cell))
        num_settlers += int(np.sum(connectivity_histogram_array))

        for ii in pdfs_to_plot_indices:
            pdf = np.copy(connectivity_histogram_array[ii])
           
            connectivity_pdf_list[ii] += pdf
            release_counts_per_cell_list[ii] += release_counts_per_cell

    else:
        continue




for pdf_index in pdfs_to_plot_indices:

#    print(pdf_index)

    #print(f"{np.shape(connectivity_pdf_list[pdf_index])}")
    #print(f"{np.shape(release_counts_per_cell_list[pdf_index])}")
    #print(f"{np.shape(release_counts_per_cell_list[pdf_index].T)}")
    #print(f"{np.shape(release_counts_per_cell_list[pdf_index][:, np.newaxis])}")
    #print(f"{np.shape(release_counts_per_cell_list[pdf_index].T[:, np.newaxis])}")
    
    # recall paul's lesson; np.newaxis really means None, and this is about the pointwise expansion business... very convenient
    #####connectivity_pdf_list[pdf_index] = connectivity_pdf_list[pdf_index] / release_counts_per_cell_list[pdf_index][:, np.newaxis]
    connectivity_pdf_list[pdf_index] = connectivity_pdf_list[pdf_index] / release_counts_per_cell_list[pdf_index][0][:, np.newaxis] # Now need to access element 0
   
    # Same result, just checking
    #for ii in range(np.shape(connectivity_pdf_list[pdf_index][0])[0]):
    #    connectivity_pdf_list[pdf_index][ii,:] /= release_counts_per_cell_list[pdf_index][0][ii]

    print(f"{np.shape(connectivity_pdf_list[pdf_index])}")
    
    if np.amax(connectivity_pdf_list[pdf_index]) > pdf_max_val:
        pdf_max_val = np.amax(connectivity_pdf_list[pdf_index])
    if np.amin(np.ma.masked_invalid(connectivity_pdf_list[pdf_index])) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(connectivity_pdf_list[pdf_index]))




settle_strength = num_settlers/num_floats




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

for pdf_index in pdfs_to_plot_indices:
   
    if pdf_index > 0:
        axis_index = pdf_index - 1
    else:
        axis_index = 0

    pdf_plot = connectivity_pdf_list[pdf_index]
    #pdf_plot = connectivity_pdf_list[pdf_index][0,:,:]   # Why is this neccesary now???

    print(f"{np.shape(pdf_plot)}")


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
        conn_csv_file_name = output_csv_dir + "/" + connectivity_file_stem + f"_annual.csv"
        #conn_csv_file_name = output_csv_dir_stem + "/" + connectivity_file_stem + f"_annual.csv"
    else:
        conn_csv_file_name = output_csv_dir + "/" + connetivity_file_stem + f"_season_{pdf_index + 1}/4.csv"
        #conn_csv_file_name = output_csv_dir_stem + "/" + connetivity_file_stem + f"_season_{pdf_index + 1}/4.csv"

    ###np.savetxt('connectivity.csv', pdf_plot, delimiter = ',')
    #np.savetxt(conn_csv_file_name, pdf_plot, delimiter = ',')

    ###pcs.append(ax.pcolormesh(X,Y,pdf_plot,cmap='jet',norm=mpl.colors.Normalize(vmin=pdf_min_val,vmax=pdf_max_val)))
    #pcs.append(ax.pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=mpl.colors.Normalize(vmin=pdf_min_val,vmax=pdf_max_val)))
    #pcs.append(ax.pcolormesh(X,Y,pdf_plot.T,cmap='jet'))
    pcs.append(ax.pcolormesh(X,Y,pdf_plot.T,cmap='jet',norm=colors.LogNorm(vmin=0.001,vmax=.1),shading='auto'))
    
    ###ax.plot([axis_index,np.shape(pdf_plot)[1]-1],[0,np.shape(pdf_plot)[0]-1],color="white",linewidth=0.5)
    ax.plot([axis_index,np.shape(pdf_plot.T)[1]-1],[0,np.shape(pdf_plot.T)[0]-1],color="white",linewidth=0.5)

    ax.set_title(subplot_title)
    #ax.title.set(subplot_title)

    #if annualOnlySwitch:
    #    break




cbar_label = "probability"

cbar = plt.colorbar(pcs[0],label=cbar_label,extend="both")
#cbar = plt.colorbar(pcs[0], ax=axs.ravel(),label=cbar_label,extend="both")

ax.set_aspect('equal')

#fig.suptitle(figure_title)

#plt.figtext(0.5, 0.02, figure_title_base, ha="center")
plt.figtext(0.5, 0.02, figure_title, ha="center")

plt.savefig(figure_file)


#plt.show()



