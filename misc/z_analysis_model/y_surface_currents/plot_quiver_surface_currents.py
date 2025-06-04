
# HACK HARDCODED BS, should be saved with means, "mask" in those files is currently just the surface temp, which is wrong
# ------------------------------------------------------------------------------------------------------------
input_model_grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_noModification.npz"
#input_model_grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"
# ------------------------------------------------------------------------------------------------------------




import os
import sys
import numpy as np
import netCDF4 
import matplotlib.pyplot as plt
import argparse
from pathlib import Path


# input arguments
# --------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("surfacecurrentfile", type=str)
parser.add_argument("oneDimensionallatlonswitch", type=int)
parser.add_argument("annual0seasonal1monthly2", type=int)
parser.add_argument("seasondex", type=int)
parser.add_argument("scbightswitch", type=int)
args = parser.parse_args()

surface_current_file = os.path.realpath(args.surfacecurrentfile)
oneDimensional_lat_lon_switch = bool(args.oneDimensionallatlonswitch)
mean_type_switch = args.annual0seasonal1monthly2
season_dex = args.seasondex
sc_bight_switch = bool(args.scbightswitch)

# --------------------------------------------------------------------


release_months_list_pre = np.array(["Dec-Feb", "March-May","June-Aug","Sep-Nov","Jan-Dec"])
#if annualOnlySwitch:
#    pdfs_to_plot_indices = [4]
#else:
#    pdfs_to_plot_indices = [0,1,2,3]
#
#release_months_list = release_months_list_pre[pdfs_to_plot_indices]



output_figures_dir = os.path.join(os.getcwd(),"z_figures")
os.makedirs(output_figures_dir,exist_ok=True)

image_file_descriptor__forcing = surface_current_file.split("___")[-1].split(".")[0]


d = np.load(input_model_grid_file)
mask = d["mask_rho"]

d = np.load(surface_current_file)
lon = d["lon"]
lat = d["lat"]
#mask = d["mask"]

if oneDimensional_lat_lon_switch:
    lon,lat = np.meshgrid(lon,lat)


match mean_type_switch:
    case 0:
        image_file_descriptor__time = "annual"
        #u = d['u_surf_annual'][season_dex]
        #v = d['v_surf_annual'][season_dex]
        u = d['u_surf_annual']
        v = d['v_surf_annual']
    case 1:
        #image_file_descriptor__time = "seasonal"
        image_file_descriptor__time = release_months_list_pre[season_dex]
        u = d['u_surf_seasonal'][season_dex] 
        v = d['v_surf_seasonal'][season_dex]
#        u = d['u_surf_seasonal']
#        v = d['v_surf_seasonal']
    case 2:
        image_file_descriptor__time = "monthly"
        #u= d['u_surf_monthly'][season_dex] 
        #v= d['v_surf_monthly'][season_dex] 
        u= d['u_surf_monthly']
        v= d['v_surf_monthly']
    case _:
        sys.exit('second argument is int: annual0seasonal1monthly2')


#u_surf_monthly = d['u_surf_monthly'] 
#v_surf_monthly = d['v_surf_monthly'] 
#u_surf_seasonal = d['u_surf_seasonal'] 
#v_surf_seasonal = d['v_surf_seasonal']
#u_surf_annual = d['u_surf_annual']
#v_surf_annual = d['v_surf_annual']

#step_bight = 5
#step_full = 10

step_bight = 5
step_full = 12

if sc_bight_switch:
    image_file_descriptor__location = "sc_bight"
    v_plot = v[0::step_bight,0::step_bight]
    u_plot = u[0::step_bight,0::step_bight]
#    mask_plot = mask[0::step_bight,0::step_bight]
#    if oneDimensional_lat_lon_switch:
#        lon_plot = lon[0::step_bight]
#        lat_plot = lat[0::step_bight]
#    else:
    lon_plot = lon[0::step_bight,0::step_bight]
    lat_plot = lat[0::step_bight,0::step_bight]
else:
    image_file_descriptor__location = "step_full_domain"
    v_plot = v[0::step_full,0::step_full]
    u_plot = u[0::step_full,0::step_full]
#    mask_plot = mask[0::step_full,0::step_full]
#    if oneDimensional_lat_lon_switch:
#        lon_plot = lon[0::step_full]
#        lat_plot = lat[0::step_full]
#    else:
    lon_plot = lon[0::step_full,0::step_full]
    lat_plot = lat[0::step_full,0::step_full]

scale = 2
headwidth = 2
headlength = 1

dpi = 300

if mean_type_switch == 1:
    nrows = 2
    ncols = 2
    figsize = (6*ncols+2,6*nrows)
    fig,axs = plt.subplots(nrows=nrows,ncols=ncols, squeeze=False, figsize=figsize, gridspec_kw={'hspace':0.0001, 'wspace':0.15}, dpi=dpi)

    for season_index in [0,1,2,3]:
        axs[0,1].pcolormesh(lon,lat,mask)
        axs[0,1].quiver(lon_plot,lat_plot,u_plot,v_plot, scale=scale, headwidth=headwidth, headlength=headlength)


else:
    nrows = 1
    ncols = 1
    figsize = (6*ncols+2,6*nrows+1)
    fig,ax = plt.subplots()

    #c = ax.pcolormesh(lon_plot,lat_plot,mask_plot)
    c = ax.pcolormesh(lon,lat,mask)
    #q = ax.quiver(lon_plot,lat_plot,u_plot,v_plot, scale=6., headwidth=4, headlength=5)
    #q = ax.quiver(lon_plot,lat_plot,u_plot,v_plot, scale=10, headwidth=1, headlength=4)
    q = ax.quiver(lon_plot,lat_plot,u_plot,v_plot, scale=scale, headwidth=headwidth, headlength=headlength)

    title = f'{image_file_descriptor__forcing}_{image_file_descriptor__time}'

    if sc_bight_switch:
        title = f'{title}\nSC Bight Detail'



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




plt.suptitle(title)


save_file = os.path.join(output_figures_dir,f"surface_velocities_quiver__{image_file_descriptor__location}__{image_file_descriptor__time}.png")

#save_file = (os.path.join(output_figures_dir,
#    f"surface_velocities_quiver__{image_file_descriptor__forcing}__{image_file_descriptor__location}__{image_file_descriptor__time}.png"))

plt.savefig(save_file)
#plt.show()

