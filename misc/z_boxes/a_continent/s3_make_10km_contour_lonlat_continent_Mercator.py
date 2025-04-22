# For some islands, contours are segmented where they intersect other islands
# So, combine them...

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from skimage import measure
import scipy.interpolate as spint
from pathlib import Path


#---------------------------------------------------------------------
#---------------------------------------------------------------------
#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
cwd = os.getcwd()
output_dir = os.path.join(str(cwd),"z_output")
Path(output_dir).mkdir(parents=True, exist_ok=True)

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)

lon_field = d["lon_rho"]
lat_field = d["lat_rho"]


#---------------------------------------------------------------------
output_file = os.path.join(output_dir,'isodistance_to_coastline_lonlat_coords_Mercator_continent')

dist_field_file = os.path.join(output_dir,'dist_2_coast_field_coastline_Mercator_continent.mat')
dist_field = scipy.io.loadmat(dist_field_file)
dist_field = dist_field['dist_field']
#---------------------------------------------------------------------



RGI = spint.RegularGridInterpolator
# create interpolator to get lat/lon at isoline points
x = np.arange(np.shape(lon_field)[0])
y = np.arange(np.shape(lon_field)[1])
rgi_lon = RGI([x,y],lon_field)
rgi_lat = RGI([x,y],lat_field)


# dist_2_coast returned values in km
# So we set "10" as our countour value

distance_from_coast = 10
contours = measure.find_contours(dist_field, distance_from_coast)
#contours = measure.find_contours(np.transpose(dist_field), distance_from_coast)


# Assume all contours returned are relevant, combine them into single contour
isoline_ij = np.vstack(contours)


#SUPER HACKY - FOR SOME REASON THERE ARE THREE EXTRA BUGGED POINTS AT THE END OF THE ISOLINE, AND THEY ARE INLAND NEAR THE SF BAY (WHICH I MASKED OUT, SO THE BUG MAY BE RELATED)
debugged_trim_length = 4
#debugged_trim_length = 3

isoline_lonlat = np.zeros((np.shape(isoline_ij)[0]-debugged_trim_length, np.shape(isoline_ij)[1]))
for ii in range(np.shape(isoline_ij)[0]-debugged_trim_length):
    isoline_lonlat[ii,0] = rgi_lon((isoline_ij[ii,0], isoline_ij[ii,1]))
    isoline_lonlat[ii,1] = rgi_lat((isoline_ij[ii,0], isoline_ij[ii,1]))

#isoline_lonlat = np.zeros(np.shape(isoline_ij))
#for ii in range(np.shape(isoline_ij)[0]):
#    isoline_lonlat[ii,0] = rgi_lon((isoline_ij[ii,0], isoline_ij[ii,1]))
#    isoline_lonlat[ii,1] = rgi_lat((isoline_ij[ii,0], isoline_ij[ii,1]))





d = {}
d["isoline_lonlat"] = isoline_lonlat
np.savez(output_file, **d)











