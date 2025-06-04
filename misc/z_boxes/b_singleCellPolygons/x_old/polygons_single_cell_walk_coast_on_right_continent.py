# Note: <cp> does not seem necessary...!?!?!?

import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# Hardcoded stuff just to get things done.  
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------

#landmass = "continent"
#landmass = "vancouver_island"
#landmass = "island_1"
#landmass = "island_2"

landmass_list = ["continent","vancouver_island","island_1","island_2"]


# The starting points need to be ocean points just off the coast - ie start in the first polygon/cell
if landmass == "continent":
    ii_max = 364
    ii = 145
    jj = 362
elif landmass == "vancouver_island":
    ii_max = 393
    ii = 367
    jj = 255
elif landmass == "island_1":
    ii_max = 393
    ii = 192
    jj = 309
elif landmass == "island_2":
    ii_max = 393
    ii = 179
    jj = 329
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------




# May be some bad design here - modiyfing variables defined outside of function...?
def add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index):
#def add_to_polygon_lists(ii,jj):
    polygon_lon = []
    polygon_lat = []

    polygon_lon.append(lon_psi[ii,jj])
    polygon_lon.append(lon_psi[ii,jj-1])
    polygon_lon.append(lon_psi[ii-1,jj-1])
    polygon_lon.append(lon_psi[ii-1,jj])
    polygon_lon.append(lon_psi[ii,jj])

    polygon_lat.append(lat_psi[ii,jj])
    polygon_lat.append(lat_psi[ii,jj-1])
    polygon_lat.append(lat_psi[ii-1,jj-1])
    polygon_lat.append(lat_psi[ii-1,jj])
    polygon_lat.append(lat_psi[ii,jj])

    list_of_polygon_vertex_lonlat_lists.append(np.array([polygon_lon,polygon_lat]))
    list_of_lists_polygon_iijj_pairs.append([ii,jj])
    polygon_index += 1

output_dir = 'z_output/'

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_noInland.npz"
grid_file_plot = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)
lon_rho = d["lon_rho"]
lat_rho = d["lat_rho"]
mask_rho = d["mask_rho"]

mask_psi = d["mask_psi"]
lon_psi = d["lon_psi"]
lat_psi = d["lat_psi"]

d = np.load(grid_file_plot)
mask_rho_plot = d["mask_rho"]

#mask_rho_og = d["mask_rho"]

#mask_rho = -1*(mask_rho_og - 1)

new_mask = np.zeros((np.shape(mask_rho)[0]+2,np.shape(mask_rho)[1]+2))



new_mask[1:-1,1:-1] = mask_rho

new_mask[0:-2,2:] += mask_rho
new_mask[0:-2,1:-1] += mask_rho
new_mask[0:-2,0:-2] += mask_rho

new_mask[1:-1,2:] += mask_rho
new_mask[1:-1,1:-1] += mask_rho
new_mask[1:-1,0:-2] += mask_rho

new_mask[2:,2:] += mask_rho
new_mask[2:,1:-1] += mask_rho
new_mask[2:,0:-2] += mask_rho


mask_rho_new = new_mask[1:-1,1:-1]
mask_rho_new[mask_rho == 0] = np.nan
mask_rho_new[mask_rho_new == 10] = np.nan

polygon_index = 0
list_of_polygon_vertex_lonlat_lists = []
list_of_lists_polygon_iijj_pairs = []

#add_to_polygon_lists(ii,jj)
add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)

cp = [ii,jj]
# Make fake last_point, to the south of starting point
# (may be better choice, depending on where starting point is - want the
# fake point to enforce the desired starting direction of walk/traversal)

# Last point = "lp"
#lp = [cp[0]-1,cp[1]] # for islands, starting in a NW corner
#lp = [cp[0],cp[1]+1] # for west coast continent, assume starting on "flat" N/S coastline.  Used to assume current point to left of last point, now want to start with last point below current point.
lp = [cp[0]-1,cp[1]] # for west coast continent, assume starting on "flat" N/S coastline

retro_polygon_dex = 0

#a = 0

horrible_hack_counter = 0
horrible_hack_variable = len(list_of_lists_polygon_iijj_pairs)

while True:

    #a += 1
    #print(a)
    #print(retro_polygon_dex)
    #print(len(list_of_lists_polygon_iijj_pairs))

    no_new_neighbor_cell_switch = False

    # cp left of lp
    if cp[1] < lp[1]:
        try:
            # look down
            if np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj - 1
                else:
                    ii = ii - 1
                    jj = jj
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0

            # look left
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii + 1
                    jj = jj
                else:
                    ii = ii
                    jj = jj -1
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look up
            elif np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj + 1
                else:
                    ii = ii + 1
                    jj = jj
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look right
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii - 1
                    jj = jj
                else:
                    ii = ii
                    jj = jj +1
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            else:
                no_new_neighbor_cell_switch = True
                retro_polygon_dex -= 1
                #break
        except:
            break

    # cp above lp
    elif cp[0] > lp[0]:
        try:
            # look left
            if np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii + 1
                    jj = jj
                else:
                    ii = ii
                    jj = jj -1
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look up
            elif np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj + 1
                else:
                    ii = ii + 1
                    jj = jj
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look right
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii -1
                    jj = jj
                else:
                    ii = ii
                    jj = jj +1
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look down
            elif np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj - 1
                else:
                    ii = ii - 1
                    jj = jj
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            else:
                no_new_neighbor_cell_switch = True
                retro_polygon_dex -= 1
                print("no new neighbors")
                #break
        except:
            break

    # cp right of lp
    elif cp[1] > lp[1]:
        try:
            # look up
            if np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj + 1
                else:
                    ii = ii + 1
                    jj = jj
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look right
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii - 1
                    jj = jj
                else:
                    ii = ii
                    jj = jj +1
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look down
            elif np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj - 1
                else:
                    ii = ii - 1
                    jj = jj
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look left
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii + 1
                    jj = jj
                else:
                    ii = ii
                    jj = jj -1
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            else:
                no_new_neighbor_cell_switch = True
                retro_polygon_dex -= 1
                #break
        except:
            break

    # cp below lp
    elif cp[0] < lp[0]:
        try:
            # look right
            if np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii - 1
                    jj = jj
                else:
                    ii = ii
                    jj = jj +1
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look down
            elif np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj - 1
                else:
                    ii = ii - 1
                    jj = jj
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look left
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii + 1
                    jj = jj
                else:
                    ii = ii
                    jj = jj -1
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            # look up
            elif np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                if np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj + 1
                else:
                    ii = ii + 1
                    jj = jj
                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                retro_polygon_dex = 0
            
            else:
                no_new_neighbor_cell_switch = True
                retro_polygon_dex -= 1
                #break
        except:
            break


    cp = [ii,jj]
#    print(cp)

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# Horrible hack
# Can't figure out why this algorithm gets stuck going around islands
# -------------------------------------------------------------------------------
    if horrible_hack_variable == len(list_of_lists_polygon_iijj_pairs):
        horrible_hack_counter += 1
    horrible_hack_variable = len(list_of_lists_polygon_iijj_pairs)
    if horrible_hack_counter > 8:
        break
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

    if no_new_neighbor_cell_switch:
        print('hi')

    if no_new_neighbor_cell_switch:
        if -1*retro_polygon_dex >= len(list_of_lists_polygon_iijj_pairs):
            print('Went all the way back to first cell!')
            break
        else:
            ii,jj = list_of_lists_polygon_iijj_pairs[retro_polygon_dex]

    #if ii >= ii_max or jj >= np.shape(mask_rho_new)[1]:
    if ii >= ii_max or jj >= np.shape(mask_rho_new)[1] or -1*retro_polygon_dex >= len(list_of_lists_polygon_iijj_pairs):
        break

    # Can't figure out why not working for islands
#    if len(list_of_lists_polygon_iijj_pairs) >= 26:
#        break


coastline_file_in = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_continent/z_output/coastline_coords_Mercator_continent.npz"
d = np.load(coastline_file_in)
coastline_lonlat_continent = d["coastline_lonlat"]

coastline_file_in = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/z_output/coastline_coords_Mercator_island_number_1.npz"
d = np.load(coastline_file_in)
coastline_lonlat_island_1 = d["coastline_lonlat"]

coastline_file_in = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/z_output/coastline_coords_Mercator_island_number_2.npz"
d = np.load(coastline_file_in)
coastline_lonlat_island_2 = d["coastline_lonlat"]

coastline_file_in = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/z_output/coastline_coords_Mercator_island_number_3.npz"
d = np.load(coastline_file_in)
coastline_lonlat_island_3 = d["coastline_lonlat"]


fig, ax = plt.subplots()
m = ax.pcolormesh(lon_rho,lat_rho,mask_rho_plot,shading="nearest")
#m = ax.pcolormesh(lon_rho,lat_rho,mask_rho_new,shading="nearest")


for ii in range(len(list_of_polygon_vertex_lonlat_lists)):
    ax.plot(list_of_polygon_vertex_lonlat_lists[ii][0,:],list_of_polygon_vertex_lonlat_lists[ii][1,:], c='r')
    ax.annotate(ii, xy = [np.mean(list_of_polygon_vertex_lonlat_lists[ii][0,:]), np.mean(list_of_polygon_vertex_lonlat_lists[ii][1,:])], ha="center", va="center", weight="bold",c='k')

ax.axis('image')

###plt.colorbar(m)

plt.show()


