# Copied from original, on Tsunami:
# /home/blaughli/tracking_project_mesoscaleModernVersion/tracking_project/practice/bounding_boxes/create_boxes/continent/coastline_define_walk_psi_bl_lonlat_continent.pyl

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import ast
import os
from pathlib import Path

output_dir = "z_output/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

coastline_file_out = os.path.join(output_dir,"coastline_coords_Mercator_continent.npz")

#---------------------------------------------------------------------
d = np.load(grid_file)

psi = d["mask_psi"]
mask_u = d["mask_u"]
mask_v = d["mask_v"]
lon_psi = d["lon_psi"]
lat_psi = d["lat_psi"]




# size of domains i and j
n_i = np.shape(psi)[0]
n_j = np.shape(psi)[1]

# In Mercator, I think I should set a cutoff row/longitude for the northern end of the coastline
ii_final = 394 # Northern Vancouver Island
#ii_final = 360
#ii_final = 330

# begin algorithm...
# walk around clockwise, coast on the right, looking left first then clockwise

coast_lon = []
coast_lat = []
coast_ii = []
coast_jj = []

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
coast_ii.append(ii)
coast_jj.append(jj)

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
        coast_ii.append(ii)
        coast_jj.append(jj)
        lp = cp
        cp = [ii,jj]
#        print(cp)
#        print([n_i,n_j])
#        print('\n')


final_coordinates_lonlat = np.zeros((len(coast_lon),2))
final_coordinates_lonlat[:,0] = coast_lon
final_coordinates_lonlat[:,1] = coast_lat

final_coordinates_iijj = np.zeros((len(coast_ii),2))
final_coordinates_iijj[:,0] = coast_ii
final_coordinates_iijj[:,1] = coast_jj

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
d["coastline_lonlat"] = final_coordinates_lonlat
d["coastline_iijj"] = final_coordinates_iijj

np.savez(coastline_file_out, **d)





