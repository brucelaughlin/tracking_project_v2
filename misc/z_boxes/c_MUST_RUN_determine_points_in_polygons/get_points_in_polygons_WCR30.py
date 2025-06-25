# Error - think I need to use lat/lon, as in version 1 I just used i/j which
# depends on the grid type... ie it's wrong to use i/j from a polygon in psi
# to bound rho points... 

# THIS REALLY SHOULD BE GENERAL, 

import os
import netCDF4
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as plt_path
from pathlib import Path
import argparse

#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"
#polygon_csv_file_path = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/bounding_boxes_lonlat_Mercator_singleCoastalCells.txt"
#
#output_file_stem_lonlat = 'points_in_boxes_lon_lat_combined_Mercator_singleCellPolygons.p'
#output_file_stem_iijj = 'points_in_boxes_ii_jj_combined_Mercator_singleCellPolygons.p'

grid_file = "/home/blaughli/tracking_project_v2/grid_data/wcr30test1_grd.nc"
polygon_csv_file_path = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/bounding_boxes_lonlat_WCR30_singleCoastalCells.txt"

output_file_stem_lonlat = 'points_in_boxes_lon_lat_combined_WCR30_singleCellPolygons.p'
output_file_stem_iijj = 'points_in_boxes_ii_jj_combined_WCR30_singleCellPolygons.p'

output_dir = "z_output/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

output_file_lonlat = os.path.join(output_dir, output_file_stem_lonlat)
output_file_iijj = os.path.join(output_dir, output_file_stem_iijj)

#---------------------------------------------------------------------


#d = np.load(grid_file)
#lon_rho = d["lon_rho"]
#lat_rho = d["lat_rho"]
d = netCDF4.Dataset(grid_file)
lon_rho = np.array(d["lon_rho"])
lat_rho = np.array(d["lat_rho"])
d.close()


# create empty list to store all of the grid lat/lon and i/j pairs as tuples
points_lon_lat = []
points_iijj = []

n_i = np.shape(lon_rho)[0]    
n_j = np.shape(lon_rho)[1]

for ii in range(n_i):
    for jj in range(n_j):
        points_lon_lat.append((lon_rho[ii,jj],lat_rho[ii,jj]))
        points_iijj.append((ii,jj))


points_lon_lat = np.array(points_lon_lat)
points_iijj = np.array(points_iijj)



# Construct list of arrays of polygon vertices from CSV file input
polygon_number = 0
list_of_polygon_vertex_lonlat_arrays = []
with open(polygon_csv_file_path) as polygon_file:
   for line in polygon_file:
        line_items = line.rstrip().split(',')
        if line_items[0].isdigit():
            if int(line_items[0]) != polygon_number:
                if polygon_number > 0:
                    list_of_polygon_vertex_lonlat_arrays.append(current_polygon_vertices)
                polygon_number += 1
                current_polygon_vertices = np.array([float(line_items[3]), float(line_items[2])])
                continue
            current_polygon_vertices = np.vstack([current_polygon_vertices, [float(line_items[3]), float(line_items[2])]]) # note that Patrick stores lat first, then lon, so I switch these
# Must append the last polygon
list_of_polygon_vertex_lonlat_arrays.append(current_polygon_vertices)



list_of_arrays_of_points_in_polygons_lonlat = []
list_of_arrays_of_points_in_polygons_iijj = []

# each "box" is a 2 by n array, with the first column being "i" coordinates, 2nd being "j"

polygon_dex = 0

for polygon_vertex_lonlat_array in list_of_polygon_vertex_lonlat_arrays:
    polygon_dex += 1
    if polygon_vertex_lonlat_array is None:
        print('box {} has value "None" ..!?'.format(polygon_dex-1))
    if polygon_vertex_lonlat_array is not None:

        path = plt_path.Path(polygon_vertex_lonlat_array)
        #path = plt_path.Path(np.transpose(polygon_vertex_lonlat_array))  # Transpose still needed?
        
        points_inside_flags = path.contains_points(points_lon_lat) 
        
        points_inside = points_lon_lat[points_inside_flags]
        
        points_in_polygon_lon = []
        points_in_polygon_lat = []
        for point in points_inside:
            points_in_polygon_lon.append(point[0])
            points_in_polygon_lat.append(point[1])
        
        points_in_polygon_lon_lat = np.array([points_in_polygon_lon,points_in_polygon_lat])
        list_of_arrays_of_points_in_polygons_lonlat.append(points_in_polygon_lon_lat)
        
        points_inside_iijj = points_iijj[points_inside_flags]
        points_in_polygon_ii = []
        points_in_polygon_jj = []
        for point in points_inside_iijj:
            points_in_polygon_ii.append(point[0])
            points_in_polygon_jj.append(point[1])
        
        points_in_polygon_iijj = np.array([points_in_polygon_ii,points_in_polygon_jj])
        list_of_arrays_of_points_in_polygons_iijj.append(points_in_polygon_iijj)



# Store all seed points in a single list, for use in making a seed file

file = open(output_file_lonlat,'wb')
pickle.dump(list_of_arrays_of_points_in_polygons_lonlat,file)
file.close()

file = open(output_file_iijj,'wb')
pickle.dump(list_of_arrays_of_points_in_polygons_iijj,file)
file.close()




