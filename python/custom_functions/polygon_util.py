# update to avoid "blobbing" of islands

import numpy as np


def dfs(current_point, previous_point, visited_points, mask_rho_valid_cells_bool, mask_rho_bool_land_only):
    visited_points.append(current_point)
    neighboring_points = generate_neighbors(current_point, previous_point, mask_rho_valid_cells_bool, mask_rho_bool_land_only)
    for neighboring_point in neighboring_points:
        if neighboring_point not in visited_points:
            dfs(neighboring_point, current_point, visited_points, mask_rho_valid_cells_bool, mask_rho_bool_land_only)


# So, will this work?  I'm importing this file and calling its functions from another script.  "dfs" above calls "generate_neighbors"... so when that happens,
# is "generate_neighbors" in the scope of dfs? 
def generate_neighbors(current_point, previous_point, mask_rho_valid_cells_bool, mask_rho_bool_land_only):

    neighboring_points = []
   
    grid_dimensions = np.shape(mask_rho_valid_cells_bool)
    
    delta_x = current_point[1] - previous_point[1]
    delta_y = current_point[0] - previous_point[0]

    delta_x_right = delta_y
    delta_y_right = -1 * delta_x
    neighbor_right = (current_point[0] + delta_y_right, current_point[1] + delta_x_right)
    
    delta_x_forward = delta_x
    delta_y_forward = delta_y
    neighbor_forward = (current_point[0] + delta_y_forward, current_point[1] + delta_x_forward)
    
    delta_x_left = -1 * delta_x_right
    delta_y_left = -1 * delta_y_right
    neighbor_left = (current_point[0] + delta_y_left, current_point[1] + delta_x_left)
    
    delta_x_backward = -1 * delta_x_forward
    delta_y_backward = -1 * delta_y_forward
    
    # LOOK RIGHT
    neighbor_right_one_forward = (neighbor_right[0]+delta_y_forward,neighbor_right[1]+delta_x_forward)
    neighbor_right_one_backward = (neighbor_right[0]+delta_y_backward,neighbor_right[1]+delta_x_backward)
    current_point_one_forward = (current_point[0]+delta_y_forward,current_point[1]+delta_x_forward)
    current_point_one_backward = (current_point[0]+delta_y_backward,current_point[1]+delta_x_backward)

    if (point_within_domain(neighbor_right, grid_dimensions) and mask_rho_valid_cells_bool[neighbor_right] and 
            ((point_within_domain(neighbor_right_one_forward, grid_dimensions) and mask_rho_bool_land_only[neighbor_right_one_forward]) or 
                (point_within_domain(neighbor_right_one_backward, grid_dimensions) and mask_rho_bool_land_only[neighbor_right_one_backward]) or
                (point_within_domain(current_point_one_forward, grid_dimensions) and mask_rho_bool_land_only[current_point_one_forward]) or
                (point_within_domain(current_point_one_backward, grid_dimensions) and mask_rho_bool_land_only[current_point_one_backward]))):

        neighboring_points.append(neighbor_right)

    # LOOK FORWARD
    neighbor_forward_one_right = (neighbor_forward[0]+delta_y_right,neighbor_forward[1]+delta_x_right)
    neighbor_forward_one_left = (neighbor_forward[0]+delta_y_left,neighbor_forward[1]+delta_x_left)
    current_point_one_right = (current_point[0]+delta_y_right,current_point[1]+delta_x_right)
    current_point_one_left = (current_point[0]+delta_y_left,current_point[1]+delta_x_left)

    if (point_within_domain(neighbor_forward, grid_dimensions) and mask_rho_valid_cells_bool[neighbor_forward] and 
            ((point_within_domain(neighbor_forward_one_right, grid_dimensions) and mask_rho_bool_land_only[neighbor_forward_one_right]) or 
                (point_within_domain(neighbor_forward_one_left, grid_dimensions) and mask_rho_bool_land_only[neighbor_forward_one_left]) or
                (point_within_domain(current_point_one_right, grid_dimensions) and mask_rho_bool_land_only[current_point_one_right]) or
                (point_within_domain(current_point_one_left, grid_dimensions) and mask_rho_bool_land_only[current_point_one_left]))):
            
        neighboring_points.append(neighbor_forward)

    # LOOK LEFT
    neighbor_left_one_forward = (neighbor_left[0]+delta_y_forward,neighbor_left[1]+delta_x_forward)
    neighbor_left_one_backward = (neighbor_left[0]+delta_y_backward,neighbor_left[1]+delta_x_backward)
    current_point_one_forward = (current_point[0]+delta_y_forward,current_point[1]+delta_x_forward)
    current_point_one_backward = (current_point[0]+delta_y_backward,current_point[1]+delta_x_backward)

    if (point_within_domain(neighbor_left, grid_dimensions) and mask_rho_valid_cells_bool[neighbor_left] and 
            ((point_within_domain(neighbor_left_one_forward, grid_dimensions) and mask_rho_bool_land_only[neighbor_left_one_forward]) or 
                (point_within_domain(neighbor_left_one_backward, grid_dimensions) and mask_rho_bool_land_only[neighbor_left_one_backward]) or
                (point_within_domain(current_point_one_forward, grid_dimensions) and mask_rho_bool_land_only[current_point_one_forward]) or
                (point_within_domain(current_point_one_backward, grid_dimensions) and mask_rho_bool_land_only[current_point_one_backward]))):
   
        neighboring_points.append(neighbor_left)
    
    return neighboring_points


def point_within_domain(point,grid_dimensions):
    if point[0] >= 0 and point[0] < grid_dimensions[0]:
        if point[1] >= 0 and point[1] < grid_dimensions[1]:
            return True
    return False


def generate_single_cell_polygons_from_list_of_points(list_of_lists_of_polygon_vertex_coordinates,visited_points,lon_psi,lat_psi):

    for point in visited_points:
        
        polygon_lons = []
        polygon_lats = []

        ii = point[0]
        jj = point[1]

        polygon_lons.append(lon_psi[ii,jj])
        polygon_lons.append(lon_psi[ii,jj-1])
        polygon_lons.append(lon_psi[ii-1,jj-1])
        polygon_lons.append(lon_psi[ii-1,jj])
        polygon_lons.append(lon_psi[ii,jj])

        polygon_lats.append(lat_psi[ii,jj])
        polygon_lats.append(lat_psi[ii,jj-1])
        polygon_lats.append(lat_psi[ii-1,jj-1])
        polygon_lats.append(lat_psi[ii-1,jj])
        polygon_lats.append(lat_psi[ii,jj])

        list_of_lists_of_polygon_vertex_coordinates.append(np.array([polygon_lons,polygon_lats]))


def add_to_polygon_datasets(ii,jj,list_of_lists_of_polygon_vertex_coordinates,set_of_polygon_center_ij_tuples):
    polygon_lons = []
    polygon_lats = []

    polygon_lons.append(lon_psi[ii,jj])
    polygon_lons.append(lon_psi[ii,jj-1])
    polygon_lons.append(lon_psi[ii-1,jj-1])
    polygon_lons.append(lon_psi[ii-1,jj])
    polygon_lons.append(lon_psi[ii,jj])

    polygon_lats.append(lat_psi[ii,jj])
    polygon_lats.append(lat_psi[ii,jj-1])
    polygon_lats.append(lat_psi[ii-1,jj-1])
    polygon_lats.append(lat_psi[ii-1,jj])
    polygon_lats.append(lat_psi[ii,jj])

    list_of_lists_of_polygon_vertex_coordinates.append(np.array([polygon_lons,polygon_lats]))
    set_of_polygon_center_ij_tuples.add((ii,jj))


def generate_walkable_mask(mask_rho_original):
    new_mask = np.zeros((np.shape(mask_rho_original)[0]+2,np.shape(mask_rho_original)[1]+2))

    new_mask[1:-1,1:-1] = mask_rho_original

    new_mask[0:-2,2:] += mask_rho_original
    new_mask[0:-2,1:-1] += mask_rho_original
    new_mask[0:-2,0:-2] += mask_rho_original

    new_mask[1:-1,2:] += mask_rho_original
    new_mask[1:-1,1:-1] += mask_rho_original
    new_mask[1:-1,0:-2] += mask_rho_original

    new_mask[2:,2:] += mask_rho_original
    new_mask[2:,1:-1] += mask_rho_original
    new_mask[2:,0:-2] += mask_rho_original

    mask_rho_valid_cells = new_mask[1:-1,1:-1]
    
    #mask_rho_land_only = np.copy(mask_rho_valid_cells)
    
    # Bug - need to "nan out" the edges of the grid... guess this means we can't use those cells, even if 
    # they're valid parts of a coast or coastal polygon.  Oh well, we don't like things at the edges anyway...
    mask_rho_valid_cells[0,:] = np.nan
    mask_rho_valid_cells[-1,:] = np.nan
    mask_rho_valid_cells[:,0] = np.nan
    mask_rho_valid_cells[:,-1] = np.nan


    mask_rho_valid_cells[mask_rho_original == 0] = np.nan # Now unnecessary, since I'm converting to bool
    mask_rho_valid_cells[mask_rho_valid_cells == 10] = np.nan

    mask_rho_valid_cells_keep = np.copy(mask_rho_valid_cells)

    mask_rho_valid_cells[(np.isnan(mask_rho_valid_cells))] = 0
    mask_rho_valid_cells_bool = np.array(mask_rho_valid_cells, dtype=bool)
    
    #mask_rho_bool_land_only = np.array(mask_rho_original, dtype=bool)
    mask_rho_bool_land_only = np.logical_not(np.array(mask_rho_original, dtype=bool))

    #return mask_rho_valid_cells_bool, mask_rho_bool_land_only
    return mask_rho_valid_cells_bool, mask_rho_bool_land_only, mask_rho_valid_cells_keep


