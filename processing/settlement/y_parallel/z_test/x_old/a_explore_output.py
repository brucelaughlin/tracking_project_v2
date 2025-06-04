import netCDF
import numpy as np
import matplotlib.pyplot as plt



file_conn = "/home/blaughli/tracking_project_v2/t_scraps/dummy_test_dir_con/v_connHist_files/conn_hist_data_a_test_file_allPLDs_version_singleFileParallel_plds_seasons_test3_bounding_boxes_lonlat_Mercator_singleCoastalCells.npz"


d = np.load(file_conn)

# In my test file where I took data for 100 floats from a real Opendrift run and then over-wrote the lat/lon values, I didn't think about particle status.  It turns out that
# 15/100 of the floats had a status of 1 (inactive) for the whole time.   

status = dset["status"] 

rc = d["release_counts_per_polygon_array"]

sc = d["prePdf_arrays_connectivity"]

release_counts_per_polygon_array = d["release_counts_per_polygon_array"]

prePdf_arrays_connectivity = d["prePdf_arrays_connectivity"]

prePdf_arrays_connectivity_noStagnation = d["prePdf_arrays_connectivity_noStagnation"]

annualDex = np.shape(prePdf_arrays_connectivity)[0] -1



