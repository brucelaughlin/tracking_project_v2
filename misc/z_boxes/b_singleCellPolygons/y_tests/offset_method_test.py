import netCDF4
import numpy as np


output_dir = 'z_output/'

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)
lon_rho = d["lon_rho"]
lat_rho = d["lat_rho"]
mask_rho = d["mask_rho"]


mask_rho = np.zeros((3,3))
mask_rho[1,1] = 1
mask_rho[0,1] = 1
mask_rho[2,1] = 1
mask_rho[1,0] = 1
mask_rho[1,2] = 1


#mask_rho = np.ones((3,3))

#new_mask = np.zeros((np.shape(mask_rho)[0]+3,np.shape(mask_rho)[1]+3))
new_mask = np.zeros((np.shape(mask_rho)[0]+2,np.shape(mask_rho)[1]+2))

new_mask[1:-1,1:-1] = mask_rho

new_mask[0:-2,2:] += mask_rho
new_mask[0:-2,1:-1] += mask_rho
new_mask[0:-2,0:-2] += mask_rho

new_mask[1:-1,2:] += mask_rho
new_mask[1:-1,1:-1] += mask_rho
new_mask[1:-1,0:-2] += mask_rho

new_mask[2:,2:] += mask_rho
new_mask[2:,1:-1] += mask_rho
new_mask[2:,0:-2] += mask_rho





