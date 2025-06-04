# Set horizontal velocities to be constant and onshore, set w=0

# velocity:
# u, ubar, v, vbar, w

file_directory = '/home/blaughli/sample_files/test_forcing/aa_Mercator/onshore_uo'

import netCDF4
import os

vars_to_modify = ['uo']
vars_to_zero = ['vo']



for filename_pre in os.listdir(file_directory):
    
    filename = os.path.join(file_directory,filename_pre)
    dset = netCDF4.Dataset(filename, 'r+')
    for var in vars_to_modify:
        dset[var][:] = 0.1
    for var in vars_to_zero:
        dset[var][:] = 0
    dset.close()










