import numpy as np
import matplotlib.pyplot as plt

# NEED RELEASE COUNTS FOR STAGNATION IGNORING PDFS

#f="/home/blaughli/tracking_project_v2/t_scraps/test_dir_con/v_connHist_files/conn_hist_data_tracking_output_configFile_000_job_00_allPLDs_version_singleFileParallel_plds_seasons_test3_bounding_boxes_lonlat_Mercator_singleCoastalCells.npz"
f="/home/blaughli/tracking_project_v2/t_scraps/test_dir_con/v_connHist_files/conn_hist_data_a_test_file_allPLDs_version_singleFileParallel_plds_seasons_test3_bounding_boxes_lonlat_Mercator_singleCoastalCells.npz"
d=np.load(f)

rca= d['release_counts_per_polygon_array']
#plt.plot(rca[0,:])
#plt.show()

pdf = d['prePdf_arrays_connectivity']
pdf_ns = d['prePdf_arrays_connectivity_noStagnation']

n_diag = 0
n_diag_ns = 0

pld_dex = np.shape(pdf)[1] -1

for ii in range(np.shape(pdf)[2]):
    n_diag += pdf[-1,pld_dex,ii,ii]
    n_diag_ns += pdf[-1,pld_dex,ii,ii]

