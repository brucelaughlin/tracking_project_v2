
import netCDF4
import numpy as np

f = "/home/blaughli/tracking_project_v2/t_scraps/dummy_test_dir_con/a_test_file.nc"

dset = netCDF4.Dataset(f, "r")

lon = dset["lon"][:]

status = dset["status"][:]

for ulon in np.unique(lon[:,125]).data:
    print(f"lon: {ulon:>6}, count: {np.sum(lon[:,128] == ulon)}")

for t in range(np.shape(status)[1]):
    print(100-np.sum(status[:,t]))
