# nearest neighbor tests

# v3: new approach to mask... i'm dumb...

import pickle
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial


base_path = '/home/blaughli/tracking_project/'
grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
grid_path_in = base_path + grid_directory + grid_file_in

mask_output_file_dir = 'practice/onshore_swim_work/'
mask_output_file_pre = 'mask_and_coord_array.p'
mask_output_file = base_path + mask_output_file_dir + mask_output_file_pre

points_type_field = 'rho'
dset = netCDF4.Dataset(grid_path_in, 'r')
mask_rho = np.array(dset['mask_rho'])
lon_field = np.array(dset['lon_{}'.format(points_type_field)])
lat_field = np.array(dset['lat_{}'.format(points_type_field)])

dset.close

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

mask_flat = mask.flatten()
#mask_flat = mask_rho.flatten()
lon_flat = lon_field.flatten()
lat_flat = lat_field.flatten()

coord_array = np.vstack([lon_flat,lat_flat])
coord_array = coord_array.T
#coord_array = np.concatenate([lon_flat,lat_flat],axis=1)

pt = [-120.6219,34.6277]
#pt = [-119.295,34.270]

distance,index = spatial.KDTree(coord_array).query(pt)

print(mask_flat[index])

#fig,ax = plt.subplots(1,2)
#ax[0].pcolormesh(lon_field,lat_field,mask_rho)
#ax[1].pcolormesh(lon_field,lat_field,mask)
#plt.show()

file = open(mask_output_file,'wb')
pickle.dump([mask_flat,coord_array],file)
file.close()



