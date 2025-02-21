#
#save_image_name = "domain_full.png"
# -----------------------------------------------------------------------------------------
#save_plot_file = save_plot_directory + save_image_name
# -----------------------------------------------------------------------------------------

# Example problematic files:
#tracking_file = '/home/blaughli/tracking_project_v2/tracking_to_do_full_nR_20_nS_15_kick_0p1/single_model_physicsOnly_1995-2020_kick_0p1/tracking_output_configFile_014_job_19.nc'

# Example normal files:
#/home/blaughli/tracking_project_v2/tracking_to_do_full_nR_20_nS_15_kick_0p1/single_model_physicsOnly_1995-2020_kick_0p1/tracking_output_configFile_014_job_19.nc


# I've been using float index -1 as input.  And with no kicks, float -1 (last float) runs aground and is stagnant for the rest of its life, which is long

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('trackingfile',type=str)
parser.add_argument('floatdex',type=int)
args = parser.parse_args()

tracking_file = args.trackingfile
float_dex = args.floatdex

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

float_mask = lon.mask[float_dex,:]
last_dex = 0
for ii in range(len(float_mask)-1,-1,-1):
    if float_mask[ii] == False:
        last_dex = ii
        break


fig, ax = plt.subplots()
#ax.axis('image')

ax.pcolormesh(lon_grid,lat_grid,mask)

ax.plot(lon[float_dex,:],lat[float_dex,:])
ax.scatter(lon[float_dex,last_dex],lat[float_dex,last_dex], c='r')

#ax.scatter(lons_final_tstep[float_dex],lats_final_tstep[float_dex], c = 'r', s=0.1)


#plt.title(title)

save_plot_file = 'float_trajectory.png'

#plt.savefig(save_plot_file)

#plt.savefig(save_plot_file, bbox_inches='tight')

#plt.show(bbox_inches='tight')
plt.show()


