import pickle
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from skimage import measure
#from geopy import distance
#from scipy import interpolate
import scipy.interpolate as spint


base_path = '/home/blaughli/tracking_project/'
grid_directory = 'grid_data/'
#grid_file_in = 'wc15_grd_no_islands.nc'
grid_file_in = 'wc15n_grd_continent.nc'
grid_path_in = base_path + grid_directory + grid_file_in
dset = netCDF4.Dataset(grid_path_in, 'r')

points_type_field = 'rho'

lon_field = np.array(dset['lon_{}'.format(points_type_field)])
lat_field = np.array(dset['lat_{}'.format(points_type_field)])

dset.close

continent_directory = 'practice/bounding_boxes/create_boxes/continent/'

output_dir = base_path + continent_directory + 'z_output/'
dist_field_file = output_dir + 'dist_2_coast_field_rho_coastline_wc15n_continent.mat'
output_file = output_dir + 'isodistance_lonlat_coords_rho_coastline_wc15n_continent.p'


#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------



dist_field = scipy.io.loadmat(dist_field_file)
dist_field = dist_field['dist_field']


plt.pcolormesh(lon_field,lat_field,dist_field.T)

#plt.scatter(lon_field[:],lat_field[:],c='y',s=0.05)

plt.show()

