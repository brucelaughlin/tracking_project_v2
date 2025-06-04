import matplotlib.pyplot as plt
import numpy as np
import netCDF4

f = "/home/cae/fiechter/WC15_era5_glorys_hindcast/wc15_fiechter_era5_glorys_avg_1995.nc"

dset = netCDF4.Dataset(f,'r')

#%run compute_mean_currents.py /home/cae/fiechter/WC15_era5_glorys_hindcast/wc15_fiechter_era5_glorys_avg_1995.nc
a = np.array(dset['u'][5,0,:,:])
b = np.array(dset['u'][5,-1,:,:])
a[a>10] = 0
a[a<-10] = 0
b[b>10] = 0
b[b<-10] = 0
fig,axs = plt.subplots(2)
axs[0].pcolormesh(a)
axs[1].pcolormesh(b)
plt.show()
