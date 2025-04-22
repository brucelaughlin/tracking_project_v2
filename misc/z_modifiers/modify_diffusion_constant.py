# looking at ncks -m wc12_his.nc, looking for diffusion coefficients and also velocity variables.


# diffusion:
# AKv:long_name = "vertical viscosity coefficient" ;
# Akv_bak:long_name = "background vertical mixing coefficient for momentum" ;


# velocity:
# u, ubar, v, vbar, w


import netCDF4
import os
import argparse

#vars_to_zero = ['u','ubar','v','vbar','w']

#file_to_change = '/home/blaughli/tracking_project_v2/misc/dummy_files/wc15n_zeroVel.nc'

#file_directory = '/home/blaughli/tracking_project_v2/misc/dummy_files/zero_vel_10days'

parser = argparse.ArgumentParser()
parser.add_argument('filedirectory', type=str)
parser.add_argument('diffusionvalue', type=str)
args = parser.parse_args()

file_directory=args.filedirectory

diffusion_value = args.diffusionvalue
#const_value = 0
#const_value = 0.0002


#vars_to_modify = ['AKv', 'Aks', 'Akt']
#vars_to_modify = ['AKv', 'AKs', 'AKt']
vars_to_modify = ['nl_visc2', 'Akk_bak', 'Akp_bak', 'AKv', 'Akv_bak']


for filename_pre in os.listdir(file_directory):
    
    filename = os.path.join(file_directory,filename_pre)
    dset = netCDF4.Dataset(filename, 'r+')

    for var in vars_to_modify:
        dset[var][:] = diffusion_value 


dset.close()










