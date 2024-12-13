
# velocity:
# u, ubar, v, vbar, w

import os
import netCDF4

def files(path):
    for file in os.listdir(path):
        if os.path.isfile(os.path.join(path,file)):
            yield file

vars_to_zero = ['u','ubar','v','vbar','w']

files_path = '/home/blaughli/tracking_project/practice/onshore_swim_work/z_dummy_model_input/'



for file in files(files_path):

    dset = netCDF4.Dataset(files_path + file, 'r+')

    for var in vars_to_zero:
        dset[var][:] = 0

    dset.close()










