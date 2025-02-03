#
#save_image_name = "domain_full.png"
# -----------------------------------------------------------------------------------------
#save_plot_file = save_plot_directory + save_image_name
# -----------------------------------------------------------------------------------------

# Test file:

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

tracking_file = '/home/blaughli/tracking_project_v2/tracking_to_do_memoryUsageTest/single_model_physicsOnly_1995-2020_kick_0p1_memory_test/tracking_output_configFile_000_job_00.nc'

#parser = argparse.ArgumentParser()
#parser.add_argument('trackingfile',type=str)
#args = parser.parse_args()
#
#tracking_file = args.trackingfile

grid_file = '/home/blaughli/tracking_project_v2/grid_data/wc15n_grd.nc'

dset = netCDF4.Dataset(grid_file, 'r')
lon_grid = np.array(dset['lon_rho'])
lat_grid = np.array(dset['lat_rho'])
mask = np.array(dset['mask_rho'])
h = np.array(dset['h'])
dset.close()



dset = netCDF4.Dataset(tracking_file)
status = dset.variables['status'][:].mask
lon = dset.variables['lon'][:].data
lat = dset.variables['lat'][:].data
z = dset.variables['z'][:].data
dset.close()

first_deactivation_dex = 0
for ii in range(1,np.shape(status)[1]):
    if np.sum(status[:,ii]) > np.sum(status[:,ii-1]):
        first_deactivation_dex = ii
        break

plot_indices = status[:,first_deactivation_dex] == True

dummy_numbers = np.array(range(len(plot_indices)))
indices = dummy_numbers[plot_indices]

#np.where(np.logical_not(np.isnan(lon[indices[0],:])))

#lons_plot = lon[plot_indices,0:first_deactivation_dex]
#lats_plot = lat[plot_indices,0:first_deactivation_dex]


#lons_final_tstep = lon[:,-1]
#lats_final_tstep = lat[:,-1]
#status_final_tstep = status[:,-1]
#status_initial_tstep = status[:,0]

#num_active_final_timestep = np.sum(np.logical_not(status_final_tstep))
#num_active_initial_timestep = np.sum(np.logical_not(status_initial_tstep))

#title_1 = f'Number of active floats at final timestep: {num_active_final_timestep:0{len(str(num_active_initial_timestep))}}\n'
#title_2 = f'Max number of active floats expected:      {num_active_initial_timestep:0{len(str(num_active_initial_timestep))}}\n'
#title = title_1 + title_2

lon_plot = lon[indices[0],first_deactivation_dex-1]
lat_plot = lat[indices[0],first_deactivation_dex-1]
z_plot = z[indices[0],first_deactivation_dex-1]

fig, ax = plt.subplots()
#ax.axis('image')
ax.pcolormesh(lon_grid,lat_grid,mask)
ax.scatter(lon_plot,lat_plot,c = 'r', s=10)
#ax.scatter(lon[indices[0],first_deactivation_dex-1],lat[indices[0],first_deactivation_dex-1], c = 'r', s=1)


#plt.title(title)

#plt.savefig(save_plot_file)
#plt.savefig(save_plot_file, bbox_inches='tight')

#plt.show(bbox_inches='tight')

plt.show()


