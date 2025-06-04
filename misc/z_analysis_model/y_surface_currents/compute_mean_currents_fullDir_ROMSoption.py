# Adding ROMS option (need to trim fields so they align)

# Also now using a single averaged file, which is created with "/home/blaughli/tracking_project_v2/misc/z_modifiers/modify_history_compute_averages_append_fields.bash"

# Hardcoded reference years for time dimension - get this info by running "ncks -m" on a model file, looking at the metadata for the <time> or <ocean_time> variable
#reference_year_ROMS = 1900
#reference_year_ROMS = 1940
#reference_year_Mercator = 1950

# REFERENCE YEARS:
# wcr30: 1940
# mercator: 1950


# HACK HARDCODED BS, should be saved with means, "mask" in those files is currently just the surface temp, which is wrong
# ------------------------------------------------------------------------------------------------------------
#input_mercator_grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_noModification.npz"
input_mercator_grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"
# ------------------------------------------------------------------------------------------------------------


import datetime
import os
import numpy as np
import netCDF4 
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
from custom_functions.roms_at_rho_trim import roms_at_rho_trim


output_dir_pre = "z_output"
output_dir = os.path.join(os.getcwd(),output_dir_pre)

os.makedirs(output_dir, exist_ok=True)

# -------------
# Input Parser
parser = argparse.ArgumentParser()
#parser.add_argument('modelaveragefile', type=str)
parser.add_argument('modelforcingdir', type=str)
parser.add_argument('referenceyear')
parser.add_argument('nonromsswitch', nargs='?', type=str)

args = parser.parse_args()

non_ROMS_switch = args.nonromsswitch
#model_average_file = os.path.realpath(args.modelaveragefile)
model_forcing_dir = os.path.realpath(args.modelforcingdir)

reference_year = int(args.referenceyear)

# Why not make it more complicated
if non_ROMS_switch is not None: 
    ROMS_switch = False
else:
    ROMS_switch = True

model_forcing_file_list = [os.path.join(model_forcing_dir,f) for f in os.listdir(model_forcing_dir) if os.path.isfile(os.path.join(model_forcing_dir,f))]
model_forcing_file_list.sort()

output_file_leaf = "mean_surface_currents___" + Path(model_forcing_dir).stem
output_file = os.path.join(output_dir,output_file_leaf)

# ----------------------------------------------------------------------------------------------------
# Setup requires loading first file from forcing directory

#dset = netCDF4.Dataset(model_average_file,'r')
dset = netCDF4.Dataset(model_forcing_file_list[0],'r')

if ROMS_switch:
    ocean_time = np.array(dset['ocean_time'])
else:
    ocean_time = np.array(dset['time'])

#if ROMS_switch:
#    date_ref = datetime.datetime(reference_year_ROMS,1,1,0,0,0)
#else:
#    date_ref = datetime.datetime(reference_year_Mercator,1,1,0,0,0)
date_ref = datetime.datetime(reference_year,1,1,0,0,0)

if ROMS_switch:
    depth_dimension_surface_index = -1
    lon = np.array(dset['lon_rho'])
    lat = np.array(dset['lat_rho'])
    mask = np.array(dset['mask_rho'])
    lon = lon[1:-1,1:-1]
    lat = lat[1:-1,1:-1]
    mask = mask[1:-1,1:-1]
    u_template, v_template = roms_at_rho_trim(np.array(dset['u'][0,depth_dimension_surface_index,:,:]),np.array(dset['v'][0,depth_dimension_surface_index,:,:]))
else:
    depth_dimension_surface_index = 0
    d = np.load(input_mercator_grid_file)
    mask = d["mask_rho"]
    lon = np.array(dset['longitude'])
    lat = np.array(dset['latitude'])
    u_template=np.array(dset['uo'][0,depth_dimension_surface_index,:,:])
    v_template=np.array(dset['vo'][0,depth_dimension_surface_index,:,:])

dset.close()

u_dims = np.shape(u_template)
v_dims = np.shape(v_template)

u_surf_monthly = np.zeros((12,u_dims[0],u_dims[1]))
v_surf_monthly = np.zeros((12,v_dims[0],v_dims[1]))
u_surf_seasonal = np.zeros((4,u_dims[0],u_dims[1]))
v_surf_seasonal = np.zeros((4,v_dims[0],v_dims[1]))
u_surf_annual = np.zeros((u_dims[0],u_dims[1]))
v_surf_annual = np.zeros((u_dims[0],u_dims[1]))

season_month_bins = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]]

monthly_counts = [0] * 12
seasonal_counts = [0] * 4
total_count = 0


# Now loop over all files in model forcing directory
# ----------------------------------------------------------------------------------------------------
file_dex = 0

for input_file in model_forcing_file_list:

    file_dex += 1

    dset = netCDF4.Dataset(input_file,'r')

    if ROMS_switch:
        ocean_time = np.array(dset['ocean_time'])
    else:
        ocean_time = np.array(dset['time'])


    for tt in range(len(ocean_time)):
       
        print(f"file {file_dex}/{len(model_forcing_file_list)}, timestep {tt+1}/{len(ocean_time)}")

        if ROMS_switch:
            date_current = date_ref + datetime.timedelta(seconds = ocean_time[tt])
        else:
            date_current = date_ref + datetime.timedelta(hours = ocean_time[tt])
        current_month = date_current.month

        if ROMS_switch:
            u_surf, v_surf = roms_at_rho_trim(np.array(dset['u'][tt,depth_dimension_surface_index,:,:]),np.array(dset['v'][tt,depth_dimension_surface_index,:,:]))
        else:
            u_surf = np.array(dset['uo'][tt,depth_dimension_surface_index,:,:])
            v_surf = np.array(dset['vo'][tt,depth_dimension_surface_index,:,:])

        '''Assigning nan to unphysical currents (ie currents with fill values assigned).
        This is very HACKY and can probably be done better using masks, which I should 
        be able to add to Mercator data (non existent in the files I currently have).'''
        u_surf[u_surf>10] = np.nan
        v_surf[v_surf>10] = np.nan
        u_surf[u_surf<-10] = np.nan
        v_surf[v_surf<-10] = np.nan

        u_surf_monthly[current_month-1,:,:] += u_surf
        v_surf_monthly[current_month-1,:,:] += v_surf
        monthly_counts[current_month-1] += 1
        
        u_surf_annual += u_surf
        v_surf_annual += v_surf
        total_count += 1

        for season_index in range(len(season_month_bins)):
            if current_month in season_month_bins[season_index]:
                u_surf_seasonal[season_index,:,:] += u_surf
                v_surf_seasonal[season_index,:,:] += v_surf
                seasonal_counts[season_index] += 1

    dset.close()
# ----------------------------------------------------------------------------------------------------

u_surf_annual /= total_count
v_surf_annual /= total_count

for month_dex in range(len(monthly_counts)):
    u_surf_monthly[month_dex,:,:] /= monthly_counts[month_dex]
    v_surf_monthly[month_dex,:,:] /= monthly_counts[month_dex]

for season_dex in range(len(seasonal_counts)):
    u_surf_seasonal[season_dex,:,:] /= seasonal_counts[season_dex]
    v_surf_seasonal[season_dex,:,:] /= seasonal_counts[season_dex]

d = {}
d['u_surf_monthly'] = u_surf_monthly
d['v_surf_monthly'] = v_surf_monthly
d['u_surf_seasonal'] = u_surf_seasonal
d['v_surf_seasonal'] = v_surf_seasonal
d['u_surf_annual'] = u_surf_annual
d['v_surf_annual'] = v_surf_annual
d['lon'] = lon
d['lat'] = lat
d['mask'] = mask

np.savez(output_file, **d)

