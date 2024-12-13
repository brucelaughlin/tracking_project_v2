# This is good exploration, but something is fundamentally wrong. look at the plots, look at what I thought I was plotting...
# got back to what i'm saving and using for advection... what the hell...?


import matplotlib.pyplot as plt
import netCDF4
import pickle
import numpy as np
from scipy import spatial


# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------
base_path = '/home/blaughli/tracking_project/'

grid_directory = 'grid_data/'
grid_file_in = 'wc15n_grd.nc'
grid_path_in = base_path + grid_directory + grid_file_in
dset = netCDF4.Dataset(grid_path_in, 'r')

points_type_field = 'rho'
points_type_line = 'psi'
lon_field = np.array(dset['lon_{}'.format(points_type_field)])
lat_field = np.array(dset['lat_{}'.format(points_type_field)])
mask_grid = np.array(dset['mask_{}'.format(points_type_field)])
dset.close

swim_data_file = '/home/blaughli/tracking_project/practice/onshore_swim_work/z_production/swim_data.p'

file = open(swim_data_file,'rb')
mask,mask_flat,coord_array,onshore_swim_component_x,onshore_swim_component_y,neg_grad_norm_map_x,neg_grad_norm_map_y = pickle.load(file)
file.close()


infile_0 = "particle_data_0.p"
infile_1 = "particle_data_1.p"

file = open(infile_0,'rb')
lon_0,lat_0 = pickle.load(file)
file.close()

file = open(infile_1,'rb')
lon_1,lat_1 = pickle.load(file)
file.close()




model_lat_avg = 38 # CHANGE IF WE NEED TO BE MORE ACCURATE (ie make specific to each particle... ugh!)

meters_per_degree_lat = 111111
meters_per_degree_lon = meters_per_degree_lat * np.cos(np.radians(model_lat_avg))

swim_speed_max = 0.1 #10cm/s

time_step_seconds = 60*60 # I hour in seconds, right?

# ------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------




lons = lon_0.reshape(-1,1)
lats = lat_0.reshape(-1,1)
lonlatpairs = np.concatenate((lons,lats), axis=1)
distances,indices = spatial.KDTree(coord_array).query(lonlatpairs)


lonlatpairs_nn = coord_array[indices]

mask_timestep = mask_flat[indices].reshape(-1,1)
onshore_swim_component_x_timestep =onshore_swim_component_x[indices].reshape(-1,1)
onshore_swim_component_y_timestep =onshore_swim_component_y[indices].reshape(-1,1)


lon_move = lons + mask_timestep*onshore_swim_component_x_timestep*swim_speed_max/meters_per_degree_lon*time_step_seconds
lat_move = lats + mask_timestep*onshore_swim_component_y_timestep*swim_speed_max/meters_per_degree_lat*time_step_seconds

#lon_move = lons + mask_flat[indices]*(onshore_swim_component_x[indices]*swim_speed_max/meters_per_degree_lon*time_step_seconds)
#lat_move = lats + mask_flat[indices]*(onshore_swim_component_y[indices]*swim_speed_max/meters_per_degree_lat*time_step_seconds)

#lon_move = lon_move.reshape(-1,1)
#lat_move = lat_move.reshape(-1,1)
lonlatpairs_advect = np.concatenate((lon_move,lat_move), axis=1)



fig,ax = plt.subplots()
mesh1 = ax.pcolormesh(lon_field,lat_field,mask_grid,shading="nearest")
ax.axis('image')


for ii in range(np.shape(lons)[0]):
    
    ax.plot([lonlatpairs[ii,0],lonlatpairs_advect[ii,0]],[lonlatpairs[ii,1],lonlatpairs_advect[ii,1]],c='m')
    #ax.plot([lonlatpairs[ii,0],lonlatpairs_nn[ii,0]],[lonlatpairs[ii,1],lonlatpairs_nn[ii,1]],c='m')
    
    ax.plot(lonlatpairs[ii,0],lonlatpairs[ii,1],'co')
    ax.plot(lonlatpairs_advect[ii,0],lonlatpairs_advect[ii,1],'ro')
    #ax.plot(lonlatpairs_nn[ii,0],lonlatpairs_nn[ii,1],'ro')
    

plt.quiver(lon_field,lat_field,neg_grad_norm_map_x.T,neg_grad_norm_map_y.T,color='k',scale=80)


plt.show()
