# Copied from "stablePlotter_update2", modified to handle multiple input files from the parallelized connectivity histogram calculation

annualOnlySwitch = True

script_version = "stablePlotter_parallel_v1"

# Also set vmax to log(0.1), to match Patricks's 2011 figures
#pdf_max_val = 0.02
#pdf_max_val = 0.01
#pdf_max_val = 0.1
#pdf_max_val = 0.2
#---------------------------------------------------------------------


import numpy as np
import argparse
import os
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("connectivitydatadir")
args = parser.parse_args()

connectivity_data_dir = args.connectivitydatadir

pdfs_to_plot_indices = [0]  # Annual-only, for now

conn_mean_list = []
conn_std_list = []
conn_release_count_list = []

# Copying some SO code to loop over files
directory = os.fsencode(connectivity_data_dir)

file_counter = 0

for connectivity_data_file_byte in os.listdir(directory):

    print(file_counter)

    connectivity_data_file = os.fsdecode(connectivity_data_file_byte)
    connectivity_data_dir = os.fsdecode(directory)

    connectivity_data_file = os.path.join(connectivity_data_dir,connectivity_data_file)

    if connectivity_data_file.endswith(".npz"):

        file_counter += 1


#        print(f"FILE: {connectivity_data_file}")
#        print(f"file counter: {file_counter}")

        if file_counter == 1:

            # Setup stuff, now that we're combining multiple histograms 
            d = np.load(connectivity_data_file)
            pld = d['pld_days']
            release_counts_per_cell = d['release_counts_per_cell']
            connectivity_histogram_array = d['pdf_arrays_connectivity']
            num_pdfs = len(connectivity_histogram_array)
            
            release_counts_per_cell_list = []
            for ii in pdfs_to_plot_indices:
                pdf = np.copy(connectivity_histogram_array[ii])

                conn_mean_list.append(np.mean(pdf))
                conn_std_list.append(np.std(pdf))
                conn_release_count_list.append(np.sum(release_counts_per_cell))
               

        
        # Re-loading on the first loop iteration seems sloppy, but I don't see the harm 
        d = np.load(connectivity_data_file)

        release_counts_per_cell = d['release_counts_per_cell']
        connectivity_histogram_array = d['pdf_arrays_connectivity']


        #num_floats += int(np.sum(release_counts_per_cell))
        #num_settlers += int(np.sum(connectivity_histogram_array))

        for ii in pdfs_to_plot_indices:
            pdf = np.copy(connectivity_histogram_array[ii])
            conn_mean_list.append(np.mean(pdf))
            conn_std_list.append(np.std(pdf))
            conn_release_count_list.append(np.sum(release_counts_per_cell))

    else:
        continue


print("report:")
print(f'{len(np.unique(conn_mean_list.append))}/{len(conn_mean_list.append)} connectivity files have unique means')
print(f'{len(np.unique(conn_std_list.append))}/{len(conn_std_list.append)} connectivity files have unique stds')




