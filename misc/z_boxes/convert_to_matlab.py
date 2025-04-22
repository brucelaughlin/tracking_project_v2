# Need to create fields for "dist2coast.m", which I previously got from netcdf grid files

import scipy
import numpy as np

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)

lon_rho = d["lon_rho"]
lat_rho = d["lat_rho"]
mask_rho = d["mask_rho"]

d = {}
d["mask_rho"] = mask_rho
d["lon_rho"] = lon_rho
d["lat_rho"] = lat_rho

scipy.io.savemat('z_output/v_test.mat', d)
