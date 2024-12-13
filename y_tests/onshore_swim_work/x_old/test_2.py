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

# Start with continent, then modify with islands
continent_directory = 'practice/bounding_boxes/create_boxes/continent/'
output_dir = base_path + continent_directory + 'z_output/'
dist_field_file = output_dir + 'dist_2_coast_field_rho_coastline_wc15n_continent.mat'
output_file = output_dir + 'isodistance_lonlat_coords_rho_coastline_wc15n_continent.p'

# Create the "dist field" from the continent distances, then modify with the island dist fields
dist_field = scipy.io.loadmat(dist_field_file)
dist_field = dist_field['dist_field']


# Now for the islands
box_dir_general = base_path + 'practice/bounding_boxes/create_boxes/'
islands_dir_pre = 'modify_islands/'
islands_dir = box_dir_general + islands_dir_pre
output_dir = islands_dir + 'z_output/'

num_islands = 8

for island_dex in range(1,num_islands+1):

    output_file = output_dir + 'isodistance_lonlat_coords_rho_coastline_wc15n_island_number_{}.p'.format(island_dex)

    dist_field_file = output_dir + 'dist_2_coast_field_rho_coastline_wc15n_island_number_{}.mat'.format(island_dex)
    dist_field_island = scipy.io.loadmat(dist_field_file)
    dist_field_island = dist_field_island['dist_field']
    
    # Store the minimum distance to land
    dist_field = np.minimum(dist_field,dist_field_island)


# Calculate the gradient
dx,dy = np.gradient(dist_field)

r,c = np.shape(dist_field)

#r = 10
#c = 10

X,Y = np.mgrid[0:r:1,0:c:1]
#Y,X = np.mgrid[0:r,0:c]

norms = np.sqrt(np.square(dx) + np.square(dy))

dirx = np.divide(dx,norms)
diry = np.divide(dy,norms)

dirx = dirx *-1
diry = diry *-1


#fig,ax = plt.subplots()

plt.pcolormesh(lon_field,lat_field,dist_field.T,vmin=0,vmax=200)
plt.quiver(lon_field,lat_field,dirx.T,diry.T, color='r',scale=80)



#plt.pcolormesh(lon_field,lat_field,diry.T,vmin=-.2,vmax=.2,cmap='coolwarm')
#plt.pcolormesh(lon_field,lat_field,diry.T,vmin=-.2,vmax=.2,cmap='bwr')
#plt.pcolormesh(lon_field,lat_field,dirx.T)

plt.colorbar()

#plt.scatter(lon_field[:],lat_field[:],c='y',s=0.05)

plt.show()

