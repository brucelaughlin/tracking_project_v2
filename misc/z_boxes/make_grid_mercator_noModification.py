# Update to allow masking of user-specified points (ie bays, etc), stored with one i/j coordinate pair per line, ie [11,45]
# in the file:  /home/blaughli/tracking_project_v2/misc/z_boxes/extra_masked_points.txt

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import ast

# ----------------
input_file = '/data04/cedwards/forcing/mercator/reanalysis12/global-reanalysis-phy-001-030-daily_1995.nc'
# ----------------

output_dir = 'z_output/'

grid_file_out = output_dir + "mercator_diy_grid_noModification.npz"

dset = netCDF4.Dataset(input_file,'r')

lon = np.array(dset['longitude'])
lat = np.array(dset['latitude'])
u = np.array(dset['uo'][0,0,:,:])
v = np.array(dset['vo'][0,0,:,:])


# Mercator has 1D lat/lon arrays, so make grid
lon_rho,lat_rho = np.meshgrid(lon,lat)

# Use temperature as a land mask, which we'll modify
mask_rho = np.array(dset['thetao'][0,0,:,:])

dset.close()

mask_rho[mask_rho>100] = np.nan
mask_rho[mask_rho<-100] = np.nan
mask_rho /= mask_rho

mask_rho[np.isnan(mask_rho)] = 0

mask_u = (mask_rho[:,0:-1] + mask_rho[:,1:])/2
mask_u[mask_u < 1] = 0

mask_v = (mask_rho[0:-1,:] + mask_rho[1:,:])/2
mask_v[mask_v < 1] = 0

mask_psi = mask_rho[:,0:-1] + mask_rho[:,1:]
mask_psi = mask_psi[0:-1,:] + mask_psi[1:,:]
mask_psi /= 4
mask_psi[mask_psi < 1] = 0

lon_psi = lon_rho[:,0:-1] + lon_rho[:,1:]
lon_psi = lon_psi[0:-1,:] + lon_psi[1:,:]
lon_psi /= 4

lat_psi = lat_rho[:,0:-1] + lat_rho[:,1:]
lat_psi = lat_psi[0:-1,:] + lat_psi[1:,:]
lat_psi /= 4


d = {}
d["mask_rho"] = mask_rho
d["mask_psi"] = mask_psi
d["mask_u"] = mask_u
d["mask_v"] = mask_v
d["lon_psi"] = lon_psi
d["lat_psi"] = lat_psi
d["lon_rho"] = lon_rho
d["lat_rho"] = lat_rho

np.savez(grid_file_out, **d)

