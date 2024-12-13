# Try to plot trajectories over the gradient field

import pickle
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from skimage import measure
import scipy.interpolate as spint


base_path = '/home/blaughli/tracking_project/'
grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
grid_path_in = base_path + grid_directory + grid_file_in

swim_data_output_file_dir = 'practice/onshore_swim_work/z_production/'
swim_data_output_file_pre = 'swim_data_simple_v1.p'
swim_data_output_file = base_path + swim_data_output_file_dir + swim_data_output_file_pre

points_type_field = 'rho'

dset = netCDF4.Dataset(grid_path_in, 'r')

mask_rho = np.array(dset['mask_rho'])
lon_field = np.array(dset['lon_{}'.format(points_type_field)])
lat_field = np.array(dset['lat_{}'.format(points_type_field)])

dset.close



# Now for the update mask, and coordinate array


r,c = np.shape(mask_rho)

rm1 = r-1
cm1 = c-1

# Make a new mask with a "buffer region" one cell wide along all coasts
mask = mask_rho.copy()
mask = mask.astype('int')
# invert the 0's and 1's
mask = np.where((mask==0)|(mask==1),mask^1,mask)

buff_west = np.zeros(np.shape(mask))
buff_west[:,0:rm1] = mask[:,1:r]

buff_east = np.zeros(np.shape(mask))
buff_east[:,1:r] = mask[:,0:rm1]

buff_south = np.zeros(np.shape(mask))
buff_south[1:c,:] = mask[0:cm1,:]

buff_north = np.zeros(np.shape(mask))
buff_north[0:cm1,:] = mask[1:c,:]

mask = mask + buff_west + buff_east + buff_south + buff_north

mask = np.divide(mask,mask,out=np.zeros_like(mask),where=mask!=0)

# invert the 0's and 1's (back again)
mask = mask.astype('int')
mask = np.where((mask==0)|(mask==1),mask^1,mask)

# Want to avoid swimming over the boundary - so just mask out the edges...
mask[0,:] = 0
mask[-1,:] = 0
mask[:,0] = 0
mask[:,-1] = 0

mask_flat = mask.flatten()
#mask_flat = mask_rho.flatten()
lon_flat = lon_field.flatten()
lat_flat = lat_field.flatten()

coord_array = np.vstack([lon_flat,lat_flat])
coord_array = coord_array.T

#plt.pcolormesh(lon_field,lat_field,mask)
#plt.pcolormesh(lon_field,lat_field,dist_field.T,vmin=0,vmax=200)
#plt.quiver(lon_field,lat_field,onshore_swim_component_x.T,onshore_swim_component_y.T, color='r',scale=80)

#plt.colorbar()


#plt.show()

onshore_swim_component_x_map = np.divide(np.ones_like(lon_field),np.sqrt(2))
onshore_swim_component_y_map = np.divide(np.ones_like(lon_field),np.sqrt(2))



onshore_swim_component_x = onshore_swim_component_x_map.flatten()
onshore_swim_component_y = onshore_swim_component_y_map.flatten()

file = open(swim_data_output_file,'wb')
pickle.dump([mask,mask_flat,coord_array,onshore_swim_component_x,onshore_swim_component_y,onshore_swim_component_x_map,onshore_swim_component_y_map],file)
file.close()



