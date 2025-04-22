
split_index = 1000

import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import argparse

grid_file = '/home/blaughli/tracking_project_v2/grid_data/wc15n_grd.nc'
dset = netCDF4.Dataset(grid_file, 'r')
lon_grid = np.array(dset['lon_rho'])
lat_grid = np.array(dset['lat_rho'])
mask = np.array(dset['mask_rho'])
h = np.array(dset['h'])
dset.close()


parser = argparse.ArgumentParser()
parser.add_argument('inputfile')
args = parser.parse_args()

input_file = args.inputfile



dset = netCDF4.Dataset(input_file, 'r')
lon = np.array(dset['lon'])
lat = np.array(dset['lat'])
z = np.array(dset['z'])
dset.close()

lon_initial_shallow = lon[0:split_index,0]
lat_initial_shallow = lat[0:split_index,0]
lon_final_shallow = lon[0:split_index,-2]
lat_final_shallow = lat[0:split_index,-2]

lon_initial_deep = lon[split_index:,0]
lat_initial_deep = lat[split_index:,0]
lon_final_deep = lon[split_index:,-2]
lat_final_deep = lat[split_index:,-2]


fig, ax = plt.subplots()
ax.pcolormesh(lon_grid,lat_grid,mask)

ax.scatter(lon_initial_shallow, lat_initial_shallow, c='g', marker='.')
ax.scatter(lon_initial_deep, lat_initial_deep, c='r', marker='.')

ax.scatter(lon_final_shallow, lat_final_shallow, c='g', marker='x')
ax.scatter(lon_final_deep, lat_final_deep, c='r', marker='x')

ax.set_xlim([min(lon[:,0]) -1, max(lon[:,0]) + 1])
ax.set_ylim([min(lat[:,0]) -1, max(lat[:,0]) + 1])

plt.show()
