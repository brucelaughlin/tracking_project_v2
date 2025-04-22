# V2: plot numbers, to use for defining which points should be cross-shore wall endpoints

import os
from pathlib import Path
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from skimage import measure
#from geopy import distance
#from scipy import interpolate
import scipy.interpolate as spint


#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)
lon_rho = d["lon_rho"]
lat_rho = d["lat_rho"]
mask_rho = d["mask_rho"]

output_dir = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/z_output"

#cwd = os.getcwd()
#output_dir = os.path.join(str(cwd),"z_output")
#Path(output_dir).mkdir(parents=True, exist_ok=True)

#---------------------------------------------------------------------
lonLat_switch = True
#lonLat_switch = False
#---------------------------------------------------------------------


#X = np.arange(-.5,int(np.shape(lon_rho)[0]-1),1)
#Y = np.arange(-.5,int(np.shape(lon_rho)[1]-1),1)
#X = np.arange(-.5,int(np.shape(lon_rho)[0]),1)
#Y = np.arange(-.5,int(np.shape(lon_rho)[1]),1)
#print(np.shape(mask_rho))
X = range(np.shape(mask_rho)[1])
Y = range(np.shape(mask_rho)[0])


fig, ax = plt.subplots()
if lonLat_switch:
    ax.pcolormesh(lon_rho,lat_rho,mask_rho,shading="nearest")
else:
    ax.pcolormesh(X,Y,mask_rho,shading="nearest")

box_dex = 0

num_islands = 3

#num_island = 6

for island_dex in range(1,num_islands+1):   
#for island_dex in range(num_island,num_island+1):   

    coastline_file_in = os.path.join(output_dir,f'coastline_coords_Mercator_island_number_{island_dex}.npz')

    d = np.load(coastline_file_in)
    coastline_lonlat = d["coastline_lonlat"]

    ax.plot(coastline_lonlat[:,0],coastline_lonlat[:,1])
    ax.scatter(coastline_lonlat[:,0],coastline_lonlat[:,1],c='b') 
    ax.scatter(coastline_lonlat[0,0],coastline_lonlat[0,1],c='r') 
    for ii in range(np.shape(coastline_lonlat)[0]):
        ax.annotate(ii, xy = [coastline_lonlat[ii,0],coastline_lonlat[ii,1]], color='white', ha="center", va="center", fontsize=15, weight="bold")
        box_dex += 1


ax.axis('image')
plt.show()







