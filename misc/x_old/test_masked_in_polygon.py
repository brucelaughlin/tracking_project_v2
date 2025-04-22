import numpy as np
#from pathlib import Path
import matplotlib.path as plt_path

box = np.array([[0,0],[1,0],[1,1],[0,1]])

a=np.full((5),np.nan)

a[0:3] = 0.5

b=np.ma.array(a,mask=np.isnan(a))

lon = b
lat = b

lon_lat = np.zeros((1,2))

lon_lat[0,0] = lon[2]
lon_lat[0,1] = lat[2]

path = plt_path.Path(box)

flag = path.contains_points(lon_lat)


lon_lat_2 = np.zeros((5,2))

lon_lat_2[:,0] = lon[:]
lon_lat_2[:,1] = lat[:]

flag_2 = path.contains_points(lon_lat_2)
