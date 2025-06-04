# just check that columns do actually add to 1

import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as tkr
import argparse
import os
from pathlib import Path

base_path = "/home/blaughli/tracking_project/practice/bounding_boxes/final_locations/z_output/z_pre_swap/z_swapped/"

pdf_raw_file = base_path + "pdf_data_seasonal_ranges_O2_pH_baseYear_1988_WC15N_1988-2010_nRunsPerNode_15_nSeed_020_dielMigration_v1_pld_20_29_swapped.npz"

d = np.load(pdf_raw_file)

pdf_arrays_connectivity = d['pdf_arrays_connectivity']
tick_positions = d['tick_positions']
tick_labels = d['tick_labels']
first_continent_box_dex = d['first_continent_box_num']

# Normalize the histograms along columns (err... rows??) to make connectivity PDFs (for inverse, normalize along rows (err... columns???))

pdf_list = []

pdf_max_val = -999999
pdf_min_val = 999999

for hist in pdf_arrays_connectivity:
    pdf = np.copy(hist)
    row_sums = pdf.sum(axis=1)
    pdf = pdf / row_sums[:, np.newaxis]
    #pdf = np.log10(pdf)
    pdf_list.append(pdf)
    if np.amax(pdf) > pdf_max_val:
        pdf_max_val = np.amax(pdf)
    if np.amin(np.ma.masked_invalid(pdf)) < pdf_min_val:
        pdf_min_val = np.amin(np.ma.masked_invalid(pdf))

##print("real data row sums:")
#for pdf in pdf_list:
#    print(np.sum(pdf,axis=1))

a=np.array([[.05,.25,.7],[.075,.175,.75],[.1,.01,.89]])
#a=np.array([[.1,.3,.6],[.1,.3,.6],[.1,.3,.6]])

print("test data:")
print(a)

print("test data row sums:")
print(np.sum(a,axis=1))

plt.pcolormesh(a.T)
plt.colorbar()
plt.show()

#/home/blaughli/tracking_project/practice/bounding_boxes/final_locations/z_output/z_pre_swap/z_swapped/pdf_data_seasonal_ranges_O2_pH_baseYear_1988_WC15N_1988-2010_nRunsPerNode_15_nSeed_020_dielMigration_v1_pld_20_29_swapped.npz

