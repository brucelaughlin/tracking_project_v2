
import os
import numpy as np
import netCDF4 
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

input_file = '/home/blaughli/symbolic_links_ROMS/WC15_hindcast_from_atlantic_202410_unpacked_singleYear/avg_1995.nc'

#grid_file = '/home/blaughli/tracking_project_v2/grid_data/wc15n_grd.nc'
grid_file = '/home/blaughli/tracking_project_v2/grid_data/wc15_grd.nc.0'

vel_max = 100



input_dir_stem = Path(os.path.dirname(os.path.realpath(input_file))).stem


#parser = argparse.ArgumentParser()
#parser.add_argument('filepath', type=str)
#args = parser.parse_args()
#input_file = args.filepath

dset = netCDF4.Dataset(grid_file,'r')
mask_rho = np.array(dset['mask_rho'])
mask_v = np.array(dset['mask_v'])
h = np.array(dset['h'])
dset.close()

masked_dex = mask_v == 0
mask_v[masked_dex] = np.nan

masked_dex = mask_rho == 0
mask_rho[masked_dex] = np.nan

h_masked = h * mask_rho

dset = netCDF4.Dataset(input_file,'r')
lon_rho = np.array(dset['lon_rho'])
lat_rho = np.array(dset['lat_rho'])
lon_v = np.array(dset['lon_v'])
lat_v = np.array(dset['lat_v'])
u=np.array(dset['u'])
v=np.array(dset['v'])
dset.close()
#vs = v[0,-1,:,:]
#vs_1 = vs[:,0:-1]
#vs_2 = vs[:,1:]
#vs_mean = .5 * (vs_1 + vs_2)
#
#us = u[0,-1,:,:]
#us_1 = us[0:-1:]
#us_2 = us[1:,:]
#us_mean = .5 * (us_1 + us_2)

vs = v[0,-1,:,:]
vs[vs>vel_max] = 0
vs_1 = vs[0:-1,:]
vs_2 = vs[1:,:]
vs_mean = .5 * (vs_1 + vs_2)

us = u[0,-1,:,:]
us[us>vel_max] = 0
us_1 = us[:,0:-1]
us_2 = us[:,1:]
us_mean = .5 * (us_1 + us_2)

vs_mean_plot = vs_mean[:,1:-1]
us_mean_plot = us_mean[1:-1,:]

lon_plot = lon_rho[1:-1,1:-1]
lat_plot = lat_rho[1:-1,1:-1]

h_plot = h_masked[1:-1,1:-1]

mask_plot = mask_rho[1:-1,1:-1]

vs_mean_plot = vs_mean_plot * mask_plot
us_mean_plot = us_mean_plot * mask_plot

#vs_mean_plot_sub = vs_mean_plot[0::5,0::5]
#us_mean_plot_sub = us_mean_plot[0::5,0::5]
#lon_plot_sub = lon_plot[0::5,0::5]
#lat_plot_sub = lat_plot[0::5,0::5]


#lon_min = -122
#lon_max = -119
#lat_min = 32
#lat_max = 36


swap_switch = False
#swap_switch = True

#sc_bight_switch = False
sc_bight_switch = True

if sc_bight_switch:
    vs_mean_plot_sub = vs_mean_plot[0::5,0::5]
    us_mean_plot_sub = us_mean_plot[0::5,0::5]
    lon_plot_sub = lon_plot[0::5,0::5]
    lat_plot_sub = lat_plot[0::5,0::5]
else:
    vs_mean_plot_sub = vs_mean_plot[0::10,0::10]
    us_mean_plot_sub = us_mean_plot[0::10,0::10]
    lon_plot_sub = lon_plot[0::10,0::10]
    lat_plot_sub = lat_plot[0::10,0::10]

fig,ax = plt.subplots()

c = ax.pcolormesh(lon_plot,lat_plot,h_plot)
if swap_switch:
    q = ax.quiver(lon_plot_sub,lat_plot_sub,vs_mean_plot_sub,us_mean_plot_sub, scale=6., headwidth = 4, headlength = 5)
    title = f'!!! U and V vector components SWAPPED !!!\n{input_dir_stem}'
else:
    q = ax.quiver(lon_plot_sub,lat_plot_sub,us_mean_plot_sub,vs_mean_plot_sub, scale=6., headwidth = 4, headlength = 5)
    title = f'{input_dir_stem}'

if sc_bight_switch:
    title = f'{title}\nSC Bight Detail'


#q = ax.quiver(lon_plot_sub,lat_plot_sub,us_mean_plot_sub,vs_mean_plot_sub, headwidth=1, scale=10, headlength=4)

ax.quiverkey(q, X=0.8, Y=.9, U=0.2, label='= 0.2m/s', labelpos='E')
#ax.quiverkey(q, X=0.45, Y=1.03, U=0.2, label='= 0.2m/s', labelpos='E')
#ax.quiverkey(q, X=0.3, Y=1.1, U=0.1, label='Quiver key, length = 0.1m/s', labelpos='E')

if sc_bight_switch:
    lon_min = -123.5
    lon_max = -118.3
    lat_min = 31.6
    lat_max = 37.2
    ax.set_xlim(lon_min,lon_max)
    ax.set_ylim(lat_min,lat_max)
#else:
#    lon_min = -125.5
#    lon_max = -120
#    lat_min = 33
#    lat_max = 43



if swap_switch:
    if sc_bight_switch:
        save_file = 'figures/swapped_sc_bight_surface_velocities.png'
    else:
        save_file = 'figures/swapped_full_surface_velocities.png'
else:
    if sc_bight_switch:
        save_file = 'figures/sc_bight_surface_velocities.png'
    else:
        save_file = 'figures/full_surface_velocities.png'


plt.suptitle(title)

#plt.show()
plt.savefig(save_file)


