import pickle
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial


base_path = '/home/blaughli/tracking_project/'
grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
grid_path_in = base_path + grid_directory + grid_file_in
dset = netCDF4.Dataset(grid_path_in, 'r')

points_type_field = 'rho'

mask_rho = np.array(dset['mask_rho'])
lon_field = np.array(dset['lon_{}'.format(points_type_field)])
lat_field = np.array(dset['lat_{}'.format(points_type_field)])

dset.close




mask_flat = mask_rho.flatten()
lon_flat = lon_field.flatten()
lat_flat = lat_field.flatten()

coord_array = np.vstack([lon_flat,lat_flat])
coord_array = coord_array.T
#coord_array = np.concatenate([lon_flat,lat_flat],axis=1)

pt = [-119.295,34.270]

distance,index = spatial.KDTree(coord_array).query(pt)







