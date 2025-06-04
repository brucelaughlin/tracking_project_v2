# For modifying a test file, after making with ncks -d trajectory,0,99,1 in.nc out.nc

import netCDF4
import os

#filename = "/home/blaughli/tracking_project_v2/t_scraps/test_dir_con/a_test_file.nc"
filename = "/home/blaughli/tracking_project_v2/t_scraps/dummy_test_dir_con/a_test_file.nc"

lons = [-123.499,-123.578,-123.582,-123.661,-123.665]
lats = [38.751,38.571,38.836,38.833,38.922]

mid_dex = 2


dset = netCDF4.Dataset(filename, 'r+')

# Resetting all status values to 0, meaning all particles are always valid.  I did not do this at first, and it led to some confusion.  Still
# concerned this might lead to oversight - my code does use status, I just hope it does so correctly
#dset["status"][:] = 0

dset["lon"][:] = lons[mid_dex]
dset["lat"][:] = lats[mid_dex]

dset["lon"][0:10,125:] = lons[0]
dset["lat"][0:10,125:] = lats[0]

dset["lon"][10:20,125:] = lons[1]
dset["lat"][10:20,125:] = lats[1]

dset["lon"][20:30,125:] = lons[3]
dset["lat"][20:30,125:] = lats[3]

dset["lon"][30:40,125:] = lons[4]
dset["lat"][30:40,125:] = lats[4]


dset.close()










