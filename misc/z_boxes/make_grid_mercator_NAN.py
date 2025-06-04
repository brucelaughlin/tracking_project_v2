# For use in making custom Opendrift landmask

import netCDF4
import numpy as np


output_dir = 'z_output/'

grid_file_out = output_dir + "mercator_diy_grid_landNAN.npz"

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"
#grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_noInland.npz"

d = np.load(grid_file)
lon_rho = d["lon_rho"]
lat_rho = d["lat_rho"]
mask_rho = d["mask_rho"]

mask_rho[mask_rho == 0] = np.nan
mask_rho[mask_rho == 1] = 0

d = {}
d["mask_rho"] = mask_rho
d["lon_rho"] = lon_rho
d["lat_rho"] = lat_rho

np.savez(grid_file_out, **d)








