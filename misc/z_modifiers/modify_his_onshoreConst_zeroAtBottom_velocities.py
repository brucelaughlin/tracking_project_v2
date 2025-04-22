# Set horizontal velocities to be constant, >0, and onshore, set w=0
# Set bottom horizontal velocities to zero

# velocity:
# u, ubar, v, vbar, w


file_directory = '/home/blaughli/sample_files/test_forcing/constant_onshore_vel_zeroAtBottom_10days'

import netCDF4
import os
import numpy as np

vars_to_modify = ['u','ubar','v','vbar']
#vars_to_zero = ['u','ubar','v','vbar','w']

constant_vel_value = 0.1

bottom_index = 0

for filename_pre in os.listdir(file_directory):
    
    filename = os.path.join(file_directory,filename_pre)
    for var in vars_to_modify:
        dset = netCDF4.Dataset(filename, 'r+')
        
        newData = np.ones(np.shape(dset[var])) * constant_vel_value
        if len(np.shape(dset[var])) == 4:
            newData[:,bottom_index,:,:] = 0
        
        
        dset[var][:] = newData
        #dset[var][:] = constant_vel_value
        #dset[var][:,bottom_index,:,:] = 0
        dset.close()
     
        #for a,b in dset.dimensions.items():
        #    print(len(b))

        #dset[var] = newData
        #dset.close() 
        
        ##for ii in range(int(np.shape(dset[var])[0])):
        #    #dset[var][ii,:,:,:] = constant_vel_value
        #    #dset[var][ii,bottom_index,:,:] = 0
        ##    dset.close()
    
    dset = netCDF4.Dataset(filename, 'r+')
    dset['w'][:] = 0
    dset.close()










