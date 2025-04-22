# Copied from original, on Tsunami:
# /home/blaughli/tracking_project_mesoscaleModernVersion/tracking_project/practice/bounding_boxes/create_boxes/continent/coastline_define_walk_psi_bl_lonlat_continent.pyl

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import ast

# ----------------
input_file = '/data04/cedwards/forcing/mercator/reanalysis12/global-reanalysis-phy-001-030-daily_1995.nc'
# ----------------

output_dir = 'z_output/'
coastline_file_out = output_dir + 'coastline_coords_psi_file_wc15n_continent'


dset = netCDF4.Dataset(input_file,'r')

lon = np.array(dset['longitude'])
lat = np.array(dset['latitude'])
u = np.array(dset['uo'][0,0,:,:])
v = np.array(dset['vo'][0,0,:,:])


# Mercator has 1D lat/lon arrays, so make grid
lon,lat = np.meshgrid(lon,lat)

# Use temperature as a land mask, which we'll modify
mask = np.array(dset['thetao'][0,0,:,:])

dset.close()

plt.pcolormesh(mask)
#plt.pcolormesh(lon,lat,mask)


plt.show()


