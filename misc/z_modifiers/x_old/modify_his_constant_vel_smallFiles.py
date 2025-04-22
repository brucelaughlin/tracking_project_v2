
import netCDF4
import os
import h5py
import xarray as xr

os.environ['HDF5_USE_FILE_LOCKING']='FALSE'

file_directory = '/home/blaughli/forcing_samples/onshore_constant_vel_10days_1990/'
#file_to_change = '/home/blaughli/tracking_project_v2/misc/dummy_files/wc15n_zeroVel.nc'


# velocity:
# u, ubar, v, vbar, w

vel_const = 0.2

vars_to_zero = ['u','ubar','v','vbar','w']


for filename_pre in os.listdir(file_directory):
    
    filename = os.path.join(file_directory,filename_pre)
    print(filename)
    dset = netCDF4.Dataset(filename, 'r+')
    #dset = xr.open_dataset("path/to/your/dataset.nc")
    for var in vars_to_zero:
        #dset[var][:].data = vel_const
        dset[var][:] = vel_const
    dset.close()










