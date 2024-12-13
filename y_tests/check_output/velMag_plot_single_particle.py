import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('tracking_output_file', type=str)
args = parser.parse_args()

tracking_output_file = args.tracking_output_file

suspect_index = 1


dt_mins = 60
dt = dt_mins * 60

model_lat_avg = 38 # CHANGE IF WE NEED TO BE MORE ACCURATE (ie make specific to each particle... ugh!)
meters_per_degree_lat = 111111
meters_per_degree_lon = meters_per_degree_lat * np.cos(np.radians(model_lat_avg))


dset = netCDF4.Dataset(tracking_output_file, 'r')

lon_all = dset.variables['lon'][:]
lat_all = dset.variables['lat'][:]

dset.close()

lon = lon_all[0,:]
lon = lon[np.logical_not(lon.mask)].data
lat = lat_all[0,:]
lat = lat[np.logical_not(lat.mask)].data

d_x = []
d_y = []
dist = []
speed_step = []

for ii in range(len(lon)-1):
    #dx.append((lon[ii+1]-lon[ii])/meters_per_degree_lon) 
    #dy.append((lat[ii+1]-lat[ii])/meters_per_degree_lat) 

    dx = (lon[ii+1]-lon[ii]) * meters_per_degree_lon 
    dy = (lat[ii+1]-lat[ii]) * meters_per_degree_lat

    d_x.append(dx)
    d_y.append(dy)
    
    dist.append(np.sqrt(dx**2 + dy**2))
    speed_step.append(np.sqrt(dx**2 + dy**2)/dt)




fig,ax = plt.subplots(4,1)

ax[0].plot(speed_step)
ax[0].set_ylabel('speed (m/s)')
ax[0].set_xlabel('timestep')

ax[1].plot(d_x)
ax[1].set_ylabel('delta x (m)')
ax[1].set_xlabel('timestep')

ax[2].plot(d_y)
ax[2].set_ylabel('delta y (m)')
ax[2].set_xlabel('timestep')

ax[3].plot(dist)
ax[3].set_ylabel('distance (m)')
ax[3].set_xlabel('timestep')

fig.suptitle('Speed along trajectory')
#fig.suptitle('Velocities along trajectory terminating deep with SF Bay')

plt.show()

