# ----------------------------------------------
# ----------------------------------------------
# Only for testing
# ----------------------------------------------
# ----------------------------------------------

dpi = 300
max_seed_depth = 20
annualSeasonIndex = 4
tick_label_csv_file_path = "/home/blaughli/tracking_project_v2/processing/plotting/pdfs/connectivity/t_text_files/tick_labels_single_cell_Mercator.txt"

subplot_label_fontsize = 6
title_fontsize = 8
tick_label_fontsize=4.5
label_fontsize=6


#connectivity_data_dir = args.connectivitydatadir
#pDrake_switch = 

# New approach: set min val to log(0.0001), so 0.0001 is our smallest colorbar tick. indicate that in colorbar tick labels
pdf_min_val = 0
pdf_min_val_log = 0.0001

pdf_max_val = 0.01
pdf_max_val_log = 0.1


# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------
#---------------------------------------------------------------------

import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os
from pathlib import Path
from glob import glob

#parser = argparse.ArgumentParser()
#parser.add_argument("trackingdir", type=str)
#parser.add_argument("annualonlyswitch", type=int)
#parser.add_argument("ignorestagnantswitch", type=int)
##parser.add_argument("logscaleswitch", type=int)

#args = parser.parse_args()

#tracking_dir = args.trackingdir
#annualOnlySwitch = bool(args.annualonlyswitch)
#ignoreStagnantSwitch = bool(args.ignorestagnantswitch)
##logSwitch = bool(args.logscaleswitch)

tracking_dir = "/data03/blaughli/tracking_output/Mercator_coastalCells_1993_2018_kickSTD_0p0___global-reanalysis-phy-001-030-daily_1993_2018" 
annualOnlySwitch = 1
ignoreStagnantSwitch = 0 



input_connectivity_dir_stem = "v_connHist_files_test5"
#input_connectivity_dir_stem = "v_connHist_files"
output_csv_dir_stem = "w_connectivity_csv_files_test5"
#output_csv_dir_stem = "w_connectivity_csv_files"
output_figures_dir_stem = "z_figures_test5"
#output_figures_dir_stem = "z_figures"

connectivity_data_dir = os.path.join(tracking_dir,input_connectivity_dir_stem)
file_list = sorted(glob(os.path.join(connectivity_data_dir,"*.npz")))

output_figures_dir = os.path.join(connectivity_data_dir,output_figures_dir_stem)
Path(output_figures_dir).mkdir(parents=True, exist_ok=True)

output_csv_dir = os.path.join(connectivity_data_dir,output_csv_dir_stem)
Path(output_csv_dir).mkdir(parents=True, exist_ok=True)

figure_title_base = Path(tracking_dir).stem


release_months_list_pre = np.array(["Dec-Feb", "March-May","June-Aug","Sep-Nov","Jan-Dec"])
#release_months_list_pre = ["Dec-Feb", "March-May","June-Aug","Sep-Nov","Jan-Dec"]

if annualOnlySwitch:
    pdfs_to_plot_indices = [4]
else:
    pdfs_to_plot_indices = [0,1,2,3]

release_months_list = release_months_list_pre[pdfs_to_plot_indices]

histogram_distances = []

    
d = np.load(file_list[0])
pld_list = d['pld_list']

num_plds = len(pld_list)

pld_dex = num_plds-1

release_cell_index = 100
#release_cell_index = 305
#release_cell_index = 300
#release_cell_index = 387
#release_cell_index = 391


for file_dex in range(len(file_list)):
   
    print(f"file {file_dex+1}/{len(file_list)}")

    d = np.load(file_list[file_dex])


    particle_distances_per_pld_allTimesteps = d['particle_distances_per_pld_allTimesteps']
    particle_release_polygons = d['particle_release_polygons']
    #polygons_settled_per_particle = d['polygons_settled_per_particle']
   

#    if ignoreStagnantSwitch:
#        particle_distances_per_pld_settlers = d['particle_distances_per_pld_settlers_noStagnation']
#    else:
#        particle_distances_per_pld_settlers = d['particle_distances_per_pld_settlers']

    #num_polygons = int(np.max(particle_release_polygons))


    histogram_particle_indices = particle_release_polygons == release_cell_index


    #histogram_distances = particle_distances_per_pld_settlers[histogram_particle_indices]
    #histogram_distances.extend(particle_distances_per_pld_settlers[histogram_particle_indices,pld_dex])
    histogram_distances.extend(particle_distances_per_pld_allTimesteps[histogram_particle_indices,pld_dex])


bins = 100

#plt.hist(np.array(histogram_distances),bins=bins)
plt.hist(histogram_distances,bins=bins)
#plt.hist(histogram_distances)
plt.show()


