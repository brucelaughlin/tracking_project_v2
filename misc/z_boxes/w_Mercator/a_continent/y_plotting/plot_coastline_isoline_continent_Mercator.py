import os
from pathlib import Path
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from skimage import measure
import scipy.interpolate as spint
import ast


#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)

lon_field = d["lon_rho"]
lat_field = d["lat_rho"]
mask_field = d["mask_rho"]

output_dir = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_continent/z_output"



#---------------------------------------------------------------------
#---------------------------------------------------------------------

coastline_file_in = os.path.join(output_dir,'coastline_coords_Mercator_continent.npz')
d = np.load(coastline_file_in)
coastline_lonlat = d["coastline_lonlat"]
        
isoline_lonlat_file_in = os.path.join(output_dir,'isodistance_to_coastline_lonlat_coords_Mercator_continent.npz')
d = np.load(isoline_lonlat_file_in)
isoline_lonlat = d["isoline_lonlat"]


fig, ax = plt.subplots()
ax.pcolormesh(lon_field,lat_field,mask_field,shading="nearest")


ax.plot(coastline_lonlat[:,0],coastline_lonlat[:,1])
ax.plot(isoline_lonlat[:,0],isoline_lonlat[:,1])

#for ii in range(np.shape(coastline)[0]):
#    ax.annotate(ii, xy = [coastline[ii,0],coastline[ii,1]], color='white', ha="center", va="center", fontsize=15, weight="bold")
#for ii in range(np.shape(isoline_lonlat)[0]):
#    ax.annotate(ii, xy = [isoline_lonlat[ii,0],isoline_lonlat[ii,1]], color='white', ha="center", va="center", fontsize=15, weight="bold")



#ax.axis('image')
#plt.show()
ax.axis('image')
plt.show()









