import netCDF4
import numpy as np
import matplotlib.pyplot as plt

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
#grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)
lon_rho = d["lon_rho"]
lat_rho = d["lat_rho"]
mask_rho = d["mask_rho"]

mask_psi = d["mask_psi"]
lon_psi = d["lon_psi"]
lat_psi = d["lat_psi"]

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


ii_max = 394
ii_min = 145

# Hardcoding starting point(s)... Took from first version using break statement, got ii/jj of first box
ii = 145
jj = 362

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


# Define lists of steps in ii/jj directions to use in algorithm (need 4 lists, depending on where current point is in relation to last point)

ii_step_cp_left = [0,]
jj_step_cp_left = [1,0,-1,0]

ii_step_cp_above = [0,-1,0,+1]
jj_step_cp_above = [1,0,-1,0]

retro_polygon_dex = 0

while True:

    print(retro_polygon_dex)

    no_new_neighbor_cell_switch = False

    # cp left of lp
    if cp[1] < lp[1]:
        try:
            # look up
            if np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                #if [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii + 1
                    jj = jj
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            # look right
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                #if [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj +1
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0

            # look down
            elif np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                #if [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii - 1
                    jj = jj
                    #add_to_polygon_lists(ii,jj)
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0

            # look left
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                #if [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
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

    # cp above lp
    elif cp[0] > lp[0]:
        try:
            # look right
            if np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                #if [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj +1
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            # look down
            elif np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                #if [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii - 1
                    jj = jj
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            # look left
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                #if [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj -1
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            # look up
            elif np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                #if [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
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

    # cp right of lp
    elif cp[1] > lp[1]:
        try:
            # look down
            if np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                #if [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii - 1
                    jj = jj
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            # look left
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                #if [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj -1
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            # look up
            elif np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                #if [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii + 1
                    jj = jj
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            # look right
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                #if [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
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

    # cp below lp
    elif cp[0] < lp[0]:
        try:
            # look left
            if np.logical_not(np.isnan(mask_rho_new[ii,jj-1])) and [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                #if [ii,jj-1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj -1
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            # look up
            elif np.logical_not(np.isnan(mask_rho_new[ii+1,jj])) and [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                #if [ii+1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii + 1
                    jj = jj
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            # look right
            elif np.logical_not(np.isnan(mask_rho_new[ii,jj+1])) and [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                #if [ii,jj+1] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii
                    jj = jj +1
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            # look down
            elif np.logical_not(np.isnan(mask_rho_new[ii-1,jj])) and [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                #if [ii-1,jj] not in list_of_lists_polygon_iijj_pairs:
                    ii = ii - 1
                    jj = jj
                    add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
                    retro_polygon_dex = 0
            
            else:
                no_new_neighbor_cell_switch = True
                retro_polygon_dex -= 1
                #break
        except:
            break



    if no_new_neighbor_cell_switch:
        if -1*retro_polygon_dex >= len(list_of_polygon_iijj_indices):
            print('Back to first cell!')
            break
        else:
            ii,jj = list_of_polygon_iijj_indices[retro_polygon_dex]




#    # Walk poleward with coastline on right
#    # Start by looking right...?
#    for look_step in range(len(ii_step_list)):
#
#        test_ii = ii + ii_step_list[look_step]
#        test_jj = jj + jj_step_list[look_step]
#        
#        if np.logical_not(np.isnan(mask_rho_new[test_ii,test_jj])):
#            if [test_ii,test_jj] not in list_of_lists_polygon_iijj_pairs:
#                ii = test_ii
#                jj = test_jj
#                add_to_polygon_lists(ii,jj,list_of_polygon_vertex_lonlat_lists,list_of_lists_polygon_iijj_pairs,polygon_index)
#
#        break

    if ii >= ii_max or jj >= np.shape(mask_rho_new)[1]:
        break



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
m = ax.pcolormesh(lon_rho,lat_rho,mask_rho_new,shading="nearest")


for ii in range(len(list_of_polygon_vertex_lonlat_lists)):
    ax.plot(list_of_polygon_vertex_lonlat_lists[ii][0,:],list_of_polygon_vertex_lonlat_lists[ii][1,:], c='k')

###ax.plot(list_of_polygon_vertex_lonlat_lists[0][0,:],list_of_polygon_vertex_lonlat_lists[0][1,:], c='k')

ax.axis('image')

###plt.colorbar(m)

plt.show()


