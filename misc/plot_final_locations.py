#
#save_image_name = "domain_full.png"
# -----------------------------------------------------------------------------------------
#save_plot_file = save_plot_directory + save_image_name
# -----------------------------------------------------------------------------------------

# Example problematic files:
#tracking_file = '/home/blaughli/tracking_project_v2/tracking_to_do_full_nR_20_nS_15_kick_0p1/single_model_physicsOnly_1995-2020_kick_0p1/tracking_output_configFile_014_job_19.nc'

# Example normal files:
#/home/blaughli/tracking_project_v2/tracking_to_do_full_nR_20_nS_15_kick_0p1/single_model_physicsOnly_1995-2020_kick_0p1/tracking_output_configFile_014_job_19.nc


import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('trackingfile',type=str)
args = parser.parse_args()

tracking_file = args.trackingfile

grid_file = '/home/blaughli/tracking_project_v2/grid_data/wc15n_grd.nc'

dset = netCDF4.Dataset(grid_file, 'r')
lon_grid = np.array(dset['lon_rho'])
lat_grid = np.array(dset['lat_rho'])
mask = np.array(dset['mask_rho'])
h = np.array(dset['h'])
dset.close()



dset = netCDF4.Dataset(tracking_file)
status = dset.variables['status'][:]
lon = dset.variables['lon'][:]
lat = dset.variables['lat'][:]
dset.close()

lons_final_tstep = lon[:,-1]
lats_final_tstep = lat[:,-1]
status_final_tstep = status[:,-1]
status_initial_tstep = status[:,0]

num_active_final_timestep = np.sum(np.logical_not(status_final_tstep))
num_active_initial_timestep = np.sum(np.logical_not(status_initial_tstep))

title_1 = f'Number of active floats at final timestep: {num_active_final_timestep:0{len(str(num_active_initial_timestep))}}\n'
title_2 = f'Max number of active floats expected:      {num_active_initial_timestep:0{len(str(num_active_initial_timestep))}}\n'
title = title_1 + title_2

fig, ax = plt.subplots()
#ax.axis('image')

ax.pcolormesh(lon_grid,lat_grid,mask)

ax.scatter(lons_final_tstep,lats_final_tstep, c = 'r', s=0.1)


plt.title(title)

#plt.savefig(save_plot_file)
#plt.savefig(save_plot_file, bbox_inches='tight')

#plt.show(bbox_inches='tight')
plt.show()


