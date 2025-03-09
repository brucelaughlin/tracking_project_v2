# Set horizontal velocities to be constant and onshore, set w=0

# velocity:
# u, ubar, v, vbar, w

file_directory = '/home/blaughli/sample_files/dummy_files/onshore_const_vel_10days'

import netCDF4
import os

vars_to_modify = ['u','ubar','v','vbar']
#vars_to_zero = ['u','ubar','v','vbar','w']



for filename_pre in os.listdir(file_directory):
    
    filename = os.path.join(file_directory,filename_pre)
    dset = netCDF4.Dataset(filename, 'r+')
    for var in vars_to_modify:
        dset[var][:] = 0.1
    dset['w'][:] = 0
    dset.close()










