# https://stackoverflow.com/questions/20165169/change-colour-of-curve-according-to-its-y-value-in-matplotlib

# Input Files
#---------------------------------------------------------------------
swim_data_file = '/home/blaughli/tracking_project/practice/onshore_swim_work/z_production/swim_data.p'
#---------------------------------------------------------------------


import netCDF4
import pickle
import numpy as np
import matplotlib.pyplot as plt




#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
base_path = '/home/blaughli/tracking_project/'

grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
grid_path_in = base_path + grid_directory + grid_file_in
dset = netCDF4.Dataset(grid_path_in, 'r')

points_type_field = 'rho'
points_type_line = 'psi'
lon_field = np.array(dset['lon_{}'.format(points_type_field)])
lat_field = np.array(dset['lat_{}'.format(points_type_field)])
mask_rho = np.array(dset['mask_{}'.format(points_type_field)])
h = np.array(dset['h'])

dset.close


file = open(swim_data_file,'rb')
mask,mask_flat,coord_array,onshore_swim_component_x,onshore_swim_component_y,onshore_swim_component_x_map,onshore_swim_component_y_map = pickle.load(file)
file.close()

mask_rho_2 = np.ma.masked_array(mask_rho, mask_rho == 1)

fig, ax = plt.subplots()
mesh1 = ax.pcolormesh(lon_field,lat_field,mask,shading="nearest",cmap='Greens')
mesh2 = ax.pcolormesh(lon_field,lat_field,mask_rho_2,shading="nearest",cmap='Blues')
#mesh2 = ax.pcolormesh(lon_field,lat_field,mask_rho,shading="nearest")
ax.axis('image')

#plt.quiver(lon_field,lat_field,onshore_swim_component_x_map.T,onshore_swim_component_y_map.T,color='r',scale=80)




plt.show()














