# Copied from original, on Tsunami:
# /home/blaughli/tracking_project_mesoscaleModernVersion/tracking_project/practice/bounding_boxes/create_boxes/continent/coastline_define_walk_psi_bl_lonlat_continent.pyl

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import ast

# ----------------
input_file = '/data04/cedwards/forcing/mercator/reanalysis12/global-reanalysis-phy-001-030-daily_1995.nc'
# ----------------

output_dir = 'z_output/'
coastline_file_out = output_dir + 'coastline_coords_psi_file_wc15n_continent'


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

#lat/lon
#highest vertex:  (-124.649989,47.95)
#lowest vertex:   (-115.847672,30.0)
#matching ii/jj:
#highest vertex:  (257, 360)
#lowest vertex:   (363, 145)

end_coords_lon = [-115.847672,-124.649989] 
end_coords_lat = [30,47.95]

end_coords_ii = [257,363]
end_coords_jj = [360,145]

fig,ax = plt.subplots(2)

ax[0].pcolormesh(mask)
ax[0].scatter(end_coords_ii,end_coords_jj)
ax[1].pcolormesh(lon,lat,mask)
ax[1].scatter(end_coords_lon,end_coords_lat)


#plt.pcolormesh(mask)
#plt.pcolormesh(lon,lat,mask)


plt.show()


