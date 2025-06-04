# Want to compute mean currents (monthly and seasonal) for plotting/analysis

import datetime
import os
import numpy as np
import netCDF4 
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

#input_file = '/data04/cedwards/forcing/mercator/reanalysis12/global-reanalysis-phy-001-030-daily_1995.nc'
#input_file = '/home/cae/fiechter/WC15_era5_glorys_hindcast/wc15_fiechter_era5_glorys_avg_1995.nc'

# -------------
# Input Parser
parser = argparse.ArgumentParser()
parser.add_argument('filepath', type=str)
parser.add_argument('nonromsswitch', nargs='?', type=str)
args = parser.parse_args()
input_file = args.filepath
non_ROMS_switch = args.nonromsswitch
# -------------

# Why not make it more complicated
if non_ROMS_switch is not None: 
    ROMS_switch = False
else:
    ROMS_switch = True


input_dir_stem = Path(os.path.dirname(os.path.realpath(input_file))).stem

# Assumes year is last part of filename.  DOES NOT handle two years in title (ie "1995_1997")
input_dir_year = input_file.split("_")[-1].split(".")[0]

output_file = "mean_currents_" + input_dir_stem + "_" + input_dir_year

# ----------------
# Input data
dset = netCDF4.Dataset(input_file,'r')

if ROMS_switch:
    ocean_time = np.array(dset['ocean_time'])
else:
    ocean_time = np.array(dset['time'])

if ROMS_switch:
    date_ref = datetime.datetime(1900,1,1,0,0,0)
else:
    date_ref = datetime.datetime(1950,1,1,0,0,0)

if ROMS_switch:
    lon = np.array(dset['lon_rho'])
    lat = np.array(dset['lat_rho'])
    mask = np.array(dset['mask_rho'])
    u_template=np.array(dset['u'][0,0,:,:])
    v_template=np.array(dset['v'][0,0,:,:])
else:
    lon = np.array(dset['longitude'])
    lat = np.array(dset['latitude'])
    mask = np.array(dset['thetao'][0,0,:,:])
    u_template=np.array(dset['uo'][0,0,:,:])
    v_template=np.array(dset['vo'][0,0,:,:])


u_dims = np.shape(u_template)
v_dims = np.shape(v_template)

u_surf_monthly = np.zeros((12,u_dims[0],u_dims[1]))
v_surf_monthly = np.zeros((12,v_dims[0],v_dims[1]))
u_surf_seasonal = np.zeros((4,u_dims[0],u_dims[1]))
v_surf_seasonal = np.zeros((4,v_dims[0],v_dims[1]))

#ubar_monthly = np.zeros((12,u_dims[0],u_dims[1]))
#vbar_monthly = np.zeros((12,v_dims[0],v_dims[1]))
#ubar_seasonal = np.zeros((4,u_dims[0],u_dims[1]))
#vbar_seasonal = np.zeros((4,v_dims[0],v_dims[1]))

season_month_bins = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]]

monthly_counts = [0] * 12
seasonal_counts = [0] * 4

#u_surf_max = 0
#u_surf_min = 1000
#v_surf_max = 0
#v_surf_min = 1000

for tt in range(len(ocean_time)):
   
    if ROMS_switch:
        date_current = date_ref + datetime.timedelta(seconds = ocean_time[tt])
    else:
        date_current = date_ref + datetime.timedelta(hours = ocean_time[tt])
    current_month = date_current.month

    #print(date_current)
    #print(current_month)
    #print("")

    if ROMS_switch:
        u_surf = np.array(dset['u'][tt,0,:,:])
        v_surf = np.array(dset['v'][tt,0,:,:])
        #u_surf = np.array(dset['u'][tt,-1,:,:])
        #v_surf = np.array(dset['v'][tt,-1,:,:])
    else:
        u_surf = np.array(dset['uo'][tt,0,:,:])
        v_surf = np.array(dset['vo'][tt,0,:,:])

    # zero-out "masked currents"
    u_surf[u_surf>10] = np.nan
    v_surf[v_surf>10] = np.nan
    u_surf[u_surf<-10] = np.nan
    v_surf[v_surf<-10] = np.nan

#    if np.max(u_surf) > u_surf_max:
#        u_surf_max = np.max(u_surf)
#    if np.max(v_surf) > v_surf_max:
#        v_surf_max = np.max(v_surf)
#    if np.min(u_surf) < u_surf_min:
#        u_surf_min = np.min(u_surf)
#    if np.min(v_surf) < v_surf_min:
#        v_surf_min = np.min(v_surf)

    #ubar = np.array(dset['ubar'][tt,:,:])
    #vbar = np.array(dset['vbar'][tt,:,:])
    #ubar_monthly[current_month,:,:] += ubar
    #vbar_monthly[current_month,:,:] += vbar

    u_surf_monthly[current_month-1,:,:] += u_surf
    v_surf_monthly[current_month-1,:,:] += v_surf
    
    monthly_counts[current_month-1] += 1

    for season_index in range(len(season_month_bins)):
        if current_month in season_month_bins[season_index]:
            #ubar_seasonal[season_index,:,:] += ubar
            #vbar_seasonal[season_index,:,:] += vbar
            
            u_surf_seasonal[season_index,:,:] += u_surf
            v_surf_seasonal[season_index,:,:] += v_surf
            
            seasonal_counts[season_index] += 1

usm = u_surf_monthly
vsm = v_surf_monthly
uss = u_surf_seasonal
vss = v_surf_seasonal

for month_dex in range(len(monthly_counts)):
    #ubar_monthly[month_dex,:,:] /= monthly_counts[month_dex]
    #vbar_monthly[month_dex,:,:] /= monthly_counts[month_dex]
    u_surf_monthly[month_dex,:,:] /= monthly_counts[month_dex]
    v_surf_monthly[month_dex,:,:] /= monthly_counts[month_dex]

for season_dex in range(len(seasonal_counts)):
    #ubar_seasonal[season_dex,:,:] /= seasonal_counts[season_dex]
    #vbar_seasonal[season_dex,:,:] /= seasonal_counts[season_dex]
    u_surf_seasonal[season_dex,:,:] /= seasonal_counts[season_dex]
    v_surf_seasonal[season_dex,:,:] /= seasonal_counts[season_dex]

d = {}
d['u_surf_monthly'] = u_surf_monthly
d['v_surf_monthly'] = v_surf_monthly
d['u_surf_seasonal'] = u_surf_seasonal
d['v_surf_seasonal'] = v_surf_seasonal
d['lon'] = lon
d['lat'] = lat
d['mask'] = mask

np.savez(output_file, **d)

