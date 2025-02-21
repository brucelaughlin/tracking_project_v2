import matplotlib.animation as animation
import os
import numpy as np
import netCDF4 
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

# -------------
# Input Parser
#parser = argparse.ArgumentParser()
#parser.add_argument('filepath', type=str)
#args = parser.parse_args()
#input_file = args.filepath
# -------------


# ---------------------------
# constants
vel_max = 100
# ---------------------------


#input_file = '/home/blaughli/symbolic_links_ROMS/WC15_hindcast_from_atlantic_202410_unpacked_singleYear/avg_1995.nc'
input_file = '/home/cae/fiechter/WC15_era5_glorys_hindcast/wc15_fiechter_era5_glorys_avg_1995.nc'

#grid_file = '/home/blaughli/tracking_project_v2/grid_data/wc15n_grd.nc'
grid_file = '/home/blaughli/tracking_project_v2/grid_data/wc15_grd.nc.0'


# ---------------------------
# switches
swap_switch = False
#swap_switch = True

sc_bight_switch = False
#sc_bight_switch = True
# ---------------------------


input_dir_stem = Path(os.path.dirname(os.path.realpath(input_file))).stem

if swap_switch:
    title = f'!!! U and V vector components SWAPPED !!!\n{input_dir_stem}'
    if sc_bight_switch:
        save_file = 'figures/swapped_sc_bight_surface_velocities.mp4'
    else:
        save_file = 'figures/swapped_full_surface_velocities.mp4'
else:
    title = f'{input_dir_stem}'
    if sc_bight_switch:
        save_file = 'figures/sc_bight_surface_velocities.mp4'
    else:
        save_file = 'figures/full_surface_velocities.mp4'


# ----------------
# Grid data
dset = netCDF4.Dataset(grid_file,'r')
lon_rho = np.array(dset['lon_rho'])
lat_rho = np.array(dset['lat_rho'])
mask_rho = np.array(dset['mask_rho'])
mask_v = np.array(dset['mask_v'])
h = np.array(dset['h'])
dset.close()

masked_dex = mask_v == 0
mask_v[masked_dex] = np.nan

masked_dex = mask_rho == 0
mask_rho[masked_dex] = np.nan

h_masked = h * mask_rho
# ----------------



# ----------------
# Input data
dset = netCDF4.Dataset(input_file,'r')

n_timesteps = dset['u'].shape[0]


# -----------
# Animation
fig,ax = plt.subplots()
    
artists = []

for tt in range(n_timesteps):

    us=np.array(dset['u'][tt,-1,:,:])
    vs=np.array(dset['v'][tt,-1,:,:])

    #for ii in range(np.shape()

    #vs = v[0,-1,:,:]
    #vs = v[-1,:,:]
    vs[vs>vel_max] = 0
    vs_1 = vs[0:-1,:]
    vs_2 = vs[1:,:]
    vs_mean = .5 * (vs_1 + vs_2)

    #us = u[0,-1,:,:]
    #us = u[-1,:,:]
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


    # Specify sampling interval, could be automated
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
    

    container = ax.pcolormesh(lon_plot,lat_plot,h_plot)
    #container += ax.quiver(lon_plot_sub,lat_plot_sub,us_mean_plot_sub,vs_mean_plot_sub, scale=6., headwidth = 4, headlength = 5)
    #q = ax.quiver(lon_plot_sub,lat_plot_sub,us_mean_plot_sub,vs_mean_plot_sub, scale=6., headwidth = 4, headlength = 5)

    artists.append(container)

    print(tt)

#    if sc_bight_switch:
#        title = f'{title}\nSC Bight Detail'



    #ax.quiverkey(q, X=0.8, Y=.9, U=0.2, label='= 0.2m/s', labelpos='E')

#    if sc_bight_switch:
#        lon_min = -123.5
#        lon_max = -118.3
#        lat_min = 31.6
#        lat_max = 37.2
    #    ax.set_xlim(lon_min,lon_max)
    #    ax.set_ylim(lat_min,lat_max)

ani = animation.ArtistAnimation(fig=fig, artists=artists,interval=1)


#    plt.suptitle(title)

plt.show()
#plt.savefig(save_file)

ani.save(filename = save_file, writer='pillow')




