# Copied from original, on Tsunami:
# /home/blaughli/tracking_project_mesoscaleModernVersion/tracking_project/practice/bounding_boxes/create_boxes/continent/coastline_define_walk_psi_bl_lonlat_continent.pyl

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import ast
import os
from pathlib import Path

# ----------------
input_file = '/data04/cedwards/forcing/mercator/reanalysis12/global-reanalysis-phy-001-030-daily_1995.nc'
# ----------------
cwd = os.getcwd()
output_dir = os.path.join(str(cwd),"z_output")
Path(output_dir).mkdir(parents=True, exist_ok=True)

coastline_file_out = os.path.join(output_dir,'coastline_coords_Mercator_continent')


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

mask[np.isnan(mask)] = 0

mask_u = (mask[:,0:-1] + mask[:,1:])/2
mask_u[mask_u < 1] = 0

mask_v = (mask[0:-1,:] + mask[1:,:])/2
mask_v[mask_v < 1] = 0

psi = mask[:,0:-1] + mask[:,1:]
psi = psi[0:-1,:] + psi[1:,:]
psi /= 4
psi[psi < 1] = 0

lon_psi = lon[:,0:-1] + lon[:,1:]
lon_psi = lon_psi[0:-1,:] + lon_psi[1:,:]
lon_psi /= 4

lat_psi = lat[:,0:-1] + lat[:,1:]
lat_psi = lat_psi[0:-1,:] + lat_psi[1:,:]
lat_psi /= 4




# size of domains i and j
n_i = np.shape(psi)[0]
n_j = np.shape(psi)[1]

# In Mercator, I think I should set a cutoff row/longitude for the northern end of the coastline
ii_final = 360
#ii_final = 330

# begin algorithm...
# walk around clockwise, coast on the right, looking left first then clockwise

coast_lon = []
coast_lat = []

# If I got i/j of the grid dims backwards, switch the order of the starting coordinates
# Start in lower right of grid (for wc12, this is land near scb)
# This assumes we start on land.  need to adjust if that's not the case

# Now using Mercator grid - need to carefully choose starting (and perhaps ending!) points...
ii = 145
#ii = 170 
#ii = 0 
jj = 363
#jj = n_j - 1

current_psi = psi[ii,jj]

while current_psi == 0:
    jj -= 1;
    current_psi = psi[ii,jj]
 
jj = jj + 1    
 
# current point = "cp"
cp = [ii,jj]



# Make fake last_point, to the south of starting point 
# (may be better choice, depending on where starting point is - want the
# fake point to enforce the desired starting direction of walk/traversal)

# Last point = "lp"
lp = [cp[0]-1,cp[1]] # for islands, starting in a NW corner
#lp = [cp[0],cp[1]+1] # for west coast continent, starting in most SW point of land mask

coast_lon.append(lon_psi[ii,jj])
coast_lat.append(lat_psi[ii,jj])

while True:
    
    # cp left of lp
    if cp[1] < lp[1]:
        try:
            # look down
            if mask_u[ii,jj] != 1 and 0 <= jj-1 < n_j:
                ii -= 1
            # look left
            elif mask_v[ii,jj] != 1 and 0 <= jj-1 < n_j:
                jj -= 1
            # look up
            elif mask_u[ii+1,jj] != 1 and 0 <= jj-1 < n_j:
                ii += 1
            # look right
            elif mask_v[ii,jj+1] != 1 and 0 <= jj-1 < n_j:
                jj += 1
            else:
                break
        except:
            break


    # cp above lp
    elif cp[0] > lp[0]:
        try:
            # look left
            if mask_v[ii,jj] != 1 and 0 <= jj-1 < n_j:
                jj -= 1
            # look up
            elif mask_u[ii+1,jj] != 1 and 0 <= jj-1 < n_j:
                ii += 1
            # look right
            elif mask_v[ii,jj+1] != 1 and 0 <= jj-1 < n_j:
                jj += 1
            # look down
            elif mask_u[ii,jj] != 1 and 0 <= jj-1 < n_j:
                ii -= 1
            else:
                break
        except:
            break

    # cp right of lp
    elif cp[1] > lp[1]:
        try:
            # look up
            if mask_u[ii+1,jj] != 1 and 0 <= jj-1 < n_j:
                ii += 1
            # look right
            elif mask_v[ii,jj+1] != 1 and 0 <= jj-1 < n_j:
                jj += 1
            # look down
            elif mask_u[ii,jj] != 1 and 0 <= jj-1 < n_j:
                ii -= 1
            # look left
            elif mask_v[ii,jj] != 1 and 0 <= jj-1 < n_j:
                jj -= 1
            else:
                break
        except:
            break

    # cp below lp
    elif cp[0] < lp[0]:
        try:
            # look right
            if mask_v[ii,jj+1] != 1 and 0 <= jj-1 < n_j:
                jj += 1
            # look down
            elif mask_u[ii,jj] != 1 and 0 <= jj-1 < n_j:
                ii -= 1
            # look left
            elif mask_v[ii,jj] != 1 and 0 <= jj-1 < n_j:
                jj -= 1
            # look up
            elif mask_u[ii+1,jj] != 1 and 0 <= jj-1 < n_j:
                ii += 1
            else:
                break
        except:
            break

    if ii >= ii_final or jj >= n_j:
        break
    else:
        coast_lon.append(lon_psi[ii,jj])
        coast_lat.append(lat_psi[ii,jj])
        lp = cp
        cp = [ii,jj]
#        print(cp)
#        print([n_i,n_j])
#        print('\n')


final_coordinates = np.zeros((len(coast_lon),2))
final_coordinates[:,0] = coast_lon
final_coordinates[:,1] = coast_lat

lon_txt_file_path = os.path.join(output_dir,"coastline_coords_Mercator_lon_continent.txt")
lat_txt_file_path = os.path.join(output_dir,"coastline_coords_Mercator_lat_continent.txt")

# write lon/lat to separate text files
with open(lon_txt_file_path, 'w') as output:
    for coord in coast_lon:
        output.write(str(coord) + "\n")
with open(lat_txt_file_path, 'w') as output:
    for coord in coast_lat:
        output.write(str(coord) + "\n")

d = {}
d["coastline_lonlat"] = final_coordinates

np.savez(coastline_file_out, **d)





