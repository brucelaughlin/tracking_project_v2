
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os

figures_save_dir = 'z_figures'
Path(figures_save_dir).mkdir(parents=True, exist_ok=True)


mean_current_file = "/home/blaughli/tracking_project_v2/misc/z_analysis/mean_currents_WC15_era5_glorys_hindcast_1995.npz"
#mean_current_file = "/home/blaughli/tracking_project_v2/misc/z_analysis/mean_currents_reanalysis12_1995.npz"

#input_dir_stem = Path(os.path.dirname(os.path.realpath(mean_current_file))).stem
input_dir_stem = mean_current_file.split('/')[-1]

d = np.load(mean_current_file)

u_surf_monthly = d['u_surf_monthly'] 
v_surf_monthly = d['v_surf_monthly'] 
u_surf_seasonal = d['u_surf_seasonal'] 
v_surf_seasonal = d['v_surf_seasonal'] 
lon = d['lon']
lat = d['lat']
mask = d['mask']

season_dex = 0

if len(np.shape(lon)) == 1:
    # Mercator data
    lon,lat = np.meshgrid(lon,lat)

#    c = ax.pcolormesh(lon,lat,u_surf_monthly[season_dex,:,:])
#    u_plot = u_surf_monthly[season_dex,:,:]
#    v_plot = v_surf_monthly[season_dex,:,:]

    lon_plot = lon
    lat_plot = lat
    
    u_seasonal_plot = u_surf_seasonal
    v_seasonal_plot = v_surf_seasonal

    mask[mask>100] = np.nan
    mask[mask<-100] = np.nan
    mask /= mask

else:
    # ROMS data

#    for season_dex in range(np.shape(u_surf_seasonal)[0]):    

    vs = v_surf_seasonal
    vs_1 = vs[:,0:-1,:]
    vs_2 = vs[:,1:,:]
    vs_mean = .5 * (vs_1 + vs_2)

    us = u_surf_seasonal
    us_1 = us[:,:,0:-1]
    us_2 = us[:,:,1:]
    us_mean = .5 * (us_1 + us_2)

    vs_mean_plot = vs_mean[:,:,1:-1]
    us_mean_plot = us_mean[:,1:-1,:]

    lon_plot = lon[1:-1,1:-1]
    lat_plot = lat[1:-1,1:-1]

    mask_plot = mask[1:-1,1:-1]

    u_seasonal_plot = us_mean_plot
    v_seasonal_plot = vs_mean_plot
    
    mask[mask==0] = np.nan
    #mask /= mask

 #   c = ax.pcolormesh(lon_plot,lat_plot,us_mean_plot[season_dex,:,:])

fig,ax = plt.subplots(2,2)

season_dex = 0

for ii in range(2):
    for jj in range(2):


        c1 = ax[ii,jj].pcolormesh(lon,lat,mask,cmap='YlOrRd')
        #c1 = ax[ii,jj].pcolormesh(lon,lat,mask)

        q = ax[ii,jj].quiver(lon_plot,lat_plot,u_seasonal_plot[season_dex,:,:],v_seasonal_plot[season_dex,:,:], scale=6., headwidth = 4, headlength = 5)
        #q = ax.quiver(lon,lat,u_plot[season_dex,:,:],v_plot[season_dex,:,:], scale=6., headwidth = 4, headlength = 5)

        ax[ii,jj].quiverkey(q, X=0.8, Y=.9, U=0.2, label='= 0.2m/s', labelpos='E')

        lon_min = -123.5
        lon_max = -118.3
        lat_min = 31.6
        lat_max = 37.2
        ax[ii,jj].set_xlim(lon_min,lon_max)
        ax[ii,jj].set_ylim(lat_min,lat_max)

        ax[ii,jj].title.set_text(f"Season {season_dex}")
        
        season_dex += 1


        #cbar = plt.colorbar(c)

fig.suptitle(f'{input_dir_stem}')
plt.show()

