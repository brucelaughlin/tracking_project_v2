import netCDF4
import numpy as np
import matplotlib.pyplot as plt


land_type = 'continent'
#land_type = 'islands'


# ----------------
input_file = '/data04/cedwards/forcing/mercator/reanalysis12/global-reanalysis-phy-001-030-daily_1995.nc'
# ----------------

output_dir = 'z_output/'
coastline_file_in = output_dir + 'coastline_coords_psi_file_wc15n_continent.npz'


dset = netCDF4.Dataset(input_file,'r')

lon = np.array(dset['longitude'])
lat = np.array(dset['latitude'])
u = np.array(dset['uo'][0,0,:,:])
v = np.array(dset['vo'][0,0,:,:])


# Mercator has 1D lat/lon arrays, so make grid
lon,lat = np.meshgrid(lon,lat)

# Use temperature as a land mask, which we'll modify
mask = np.array(dset['thetao'][0,0,:,:])

dset.close()

mask[mask>100] = np.nan
mask[mask<-100] = np.nan
mask /= mask


d = np.load(coastline_file_in)

coastlines_lonlat = d["coastline_lonlat"]


fig, ax = plt.subplots()
ax.pcolormesh(lon,lat,mask,shading="nearest")

ax.plot(coastlines_lonlat[:,0],coastlines_lonlat[:,1])

ax.axis('image')
plt.show()









