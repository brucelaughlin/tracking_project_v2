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

model_file = "/home/blaughli/sample_files/test_forcing/aa_Mercator/onshore_uo/onshore_10days.nc"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('trackingfile',type=str)
parser.add_argument('gridfile',type=str)
parser.add_argument('floatdex',type=int)
args = parser.parse_args()

tracking_file = args.trackingfile
grid_file = args.gridfile
float_dex = args.floatdex

dset = netCDF4.Dataset(model_file,'r')
u_model = np.array(dset['uo'][0,0,:,:])

#grid_file = '/home/blaughli/tracking_project_v2/grid_data/wc15n_grd.nc'

if os.path.splitext(grid_file)[-1] == ".nc":
    dset = netCDF4.Dataset(grid_file, 'r')
    lon_grid = np.array(dset['lon_rho'])
    lat_grid = np.array(dset['lat_rho'])
    mask = np.array(dset['mask_rho'])
    dset.close()
elif os.path.splitext(grid_file)[-1] == ".npz":
    dset = np.load(grid_file, 'r')
    lon_grid = np.array(dset['lon_rho'])
    lat_grid = np.array(dset['lat_rho'])
    mask = np.array(dset['mask_rho'])


dset = netCDF4.Dataset(tracking_file)
status = dset.variables['status'][:]
lon = dset.variables['lon'][:]
lat = dset.variables['lat'][:]
dset.close()

if lon.mask.any():
#if lon.mask[0]:
    print('masked trajectories present')
    float_mask = lon.mask[float_dex,:]
    last_dex = 0
    for ii in range(len(float_mask)-1,-1,-1):
        if float_mask[ii] == False:
            last_dex = ii
            break
else:
    last_dex = len(lon[float_dex,:]) - 1


fig, ax = plt.subplots()
#ax.axis('image')

#pcm = ax.pcolormesh(lon_grid,lat_grid,u_model)
ax.pcolormesh(lon_grid,lat_grid,mask)

np.random.seed(0)

for ii in range(np.shape(lon)[0]):
    print(ii)
    r = np.random.rand()
    g = np.random.rand()
    b = np.random.rand()

    ax.plot(lon[ii,:],lat[ii,:],color = (r,g,b))
    #ax.plot(lon[float_dex,:],lat[float_dex,:])
    ax.scatter([lon[ii,0],lon[ii,-1]],[lat[ii,0],lat[ii,-1]], color=(r,g,b))
    #ax.scatter(lon[float_dex,last_dex],lat[float_dex,last_dex], color=(r,g,b))
    #ax.scatter(lon[float_dex,last_dex],lat[float_dex,last_dex], c='r')

#ax.plot(lon[float_dex,:],lat[float_dex,:])
#ax.scatter(lon[float_dex,last_dex],lat[float_dex,last_dex], c='r')

#plt.legend()

#plt.title(title)

#plt.colorbar(pcm)

plt.xlim(-123.07,-121.99)
plt.ylim(37.34,38.57)

save_plot_file = 'float_trajectory.png'

#plt.savefig(save_plot_file)

#plt.savefig(save_plot_file, bbox_inches='tight')

#plt.show(bbox_inches='tight')
plt.show()


