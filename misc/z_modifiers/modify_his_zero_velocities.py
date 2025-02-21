# looking at ncks -m wc12_his.nc, looking for diffusion coefficients and also velocity variables.


# diffusion:
# AKv:long_name = "vertical viscosity coefficient" ;
# Akv_bak:long_name = "background vertical mixing coefficient for momentum" ;


# velocity:
# u, ubar, v, vbar, w


import netCDF4
import os

vars_to_zero = ['u','ubar','v','vbar','w']

#file_to_change = '/home/blaughli/tracking_project_v2/misc/dummy_files/wc15n_zeroVel.nc'

file_directory = '/home/blaughli/tracking_project_v2/misc/dummy_files/zero_vel_10days'


for filename_pre in os.listdir(file_directory):
    
    filename = os.path.join(file_directory,filename_pre)
    dset = netCDF4.Dataset(filename, 'r+')
    for var in vars_to_zero:
        dset[var][:] = 0
    dset.close()










