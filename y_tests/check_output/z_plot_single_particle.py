import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('tracking_output_file', type=str)
args = parser.parse_args()

tracking_output_file = args.tracking_output_file

#tracking_output_file = '/data/blaughli/tracking_output/baseYear_1988/WC15N_1988-2010_nRunsPerNode_15_nSeed_020_physicsOnly/tracking_output_calcDT_060_saveDT_1440_buffer_100_nSeed_020_startNudge_000000.nc'
#tracking_output_file = '/data/blaughli/tracking_output/baseYear_1988_phantomVel/WC15N_1988-2010_nSeed_020_physicsOnly/tracking_output_calcDT_060_saveDT_0060_buffer_100_nSeed_020_startNudge_000000.nc'

#suspect_index = 100378
suspect_index = 1




dset = netCDF4.Dataset(tracking_output_file, 'r')

lon_all = dset.variables['lon'][:]
lat_all = dset.variables['lat'][:]
z_all = dset.variables['z'][:]

u_all = dset.variables['x_sea_water_velocity'][:]
v_all = dset.variables['y_sea_water_velocity'][:]
w_all = dset.variables['upward_sea_water_velocity'][:]

dset.close()

u = u_all[0,:]
u = u[np.logical_not(u.mask)].data
v = v_all[0,:]
v = v[np.logical_not(v.mask)].data
w = w_all[0,:]
w = w[np.logical_not(w.mask)].data
lon = lon_all[0,:]
lon = lon[np.logical_not(lon.mask)].data
lat = lat_all[0,:]
lat = lat[np.logical_not(lat.mask)].data
z = z_all[0,:]
z = z[np.logical_not(z.mask)].data

fig,ax = plt.subplots(1,1)

ax.plot(z)
ax.set_ylabel('z (m)')
ax.set_xlabel('timestep')

fig.suptitle('Depth along trajectory')
#fig.suptitle('Velocities along trajectory terminating deep with SF Bay')

plt.show()

