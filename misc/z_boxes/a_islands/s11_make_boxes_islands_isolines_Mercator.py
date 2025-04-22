# lat/lon area calc taken from:
# https://stackoverflow.com/questions/68118907/shapely-pyproj-find-area-in-m2-of-a-polygon-created-from-latitude-and-longi

# V2: For blob islands, no longer doing the threshold method.  Cross-shore walls re pre-determined, by "islands_blob_fixed_wall_coords_V2.txt".
# May want to do the same thing for all islands!
# ALSO:  Trying new method, where there aren't "inshore" and "offshore" boxes anymore


# Using the island coastlines, including the "blob" of islands 1-4, along with the rotated isolines (ie starting points of isolines
# correspond with starting points of coastlines).
# Idea: split each island into two halves, an upper and lower, using the "bounding points" previously determined.  Then proceed
# as if doing the box calculations for two separate coastlines and isolines


# NOTE: I stored the 2nd isoline bounding point as a negative number.  So address that accordingly.  

import pickle
import os
from pathlib import Path
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from geopy.distance import great_circle
import scipy.interpolate as spint
import ast

from pyproj import Geod
#from shapely.geometry import Polygon
from shapely.geometry import LineString, Point, Polygon


#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

cwd = os.getcwd()
output_dir = os.path.join(str(cwd),"z_output")
Path(output_dir).mkdir(parents=True, exist_ok=True)

blob_coastline_north_filepath = os.path.join(output_dir,"coastline_coords_Mercator_island_1_through_2_combined_north.npz") 
blob_coastline_south_filepath = os.path.join(output_dir,"coastline_coords_Mercator_island_1_through_2_combined_south.npz") 
blob_isoline_filepath = os.path.join(output_dir,"isodistance_lonlat_coords_Mercator_island_1_through_2_blob_rotated.npz")

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)

lon_field = d["lon_rho"]
lat_field = d["lat_rho"]
mask = d["mask_rho"]

first_non_blob_island_number = 3
last_non_blob_island_number_plus_one = 4

wall_bp_file_list = []

wall_bp_list_file_blob = os.path.join(cwd,"s10_islands_fixed_wall_coords_blob.txt")
wall_bp_file_list.append(wall_bp_list_file_blob)

wall_bp_list_template = "s10_islands_fixed_wall_coords_island_{}.txt"
for island_number in range(first_non_blob_island_number,last_non_blob_island_number_plus_one):
    wall_bp_file_list.append(os.path.join(cwd,wall_bp_list_template.format(island_number)))



num_islands = 3

num_last_blob_island = 2

bp_file_dex = 0

for island_number in range(num_last_blob_island,num_islands+1):
#for island_number in range(num_last_blob_island,num_last_blob_island+1):
#for island_number in range(6,7):

    print('\n')
    print('\n')
    print('island number: {}'.format(island_number))
    print('\n')
    print('\n')

    # Load coastline and isoline coordinates
    if island_number == num_last_blob_island:
        
        isoline_file_in = blob_isoline_filepath

        # Load the coastlines
        d = np.load(blob_coastline_north_filepath)
        coastline_lonlat_1 = d["coastline_lonlat_north"]
        d = np.load(blob_coastline_south_filepath)
        coastline_lonlat_2 = d["coastline_lonlat_south"]
    

        coastline_lonlat = np.append(coastline_lonlat_1,coastline_lonlat_2[1:,:],axis=0)
        #coastline_lonlat = np.append(coastline_lonlat_1,coastline_lonlat_2,axis=0)

    else:

        isoline_file_in = os.path.join(output_dir,'isodistance_lonlat_coords_Mercator_island_number_{}_rotated.npz'.format(island_number))
        coastline_file_in = os.path.join(output_dir,'coastline_coords_Mercator_island_number_{}.npz'.format(island_number))

        # Load the coastlines
        d = np.load(coastline_file_in)
        coastline_lonlat = d["coastline_lonlat"]

    # Load the isolines
    d = np.load(isoline_file_in)
    isoline_lonlat = d["isoline_lonlat"]


    bounding_boxes_file_out = os.path.join(output_dir,'bounding_boxes_lonlat_Mercator_island_number_{}.p'.format(island_number))
    #bounding_boxes_file_out = os.path.join(output_dir,'bounding_boxes_lonlat_Mercator_island_number_{}'.format(island_number))
    #bounding_boxes_file_out = output_dir + 'bounding_boxes_lonlat_Mercator_island_number_{}.p'.format(island_number)

    coast_lon = list(coastline_lonlat[:,0])
    coast_lat = list(coastline_lonlat[:,1])
    isoline_lon = list(isoline_lonlat[:,0])
    isoline_lat = list(isoline_lonlat[:,1])
    
    num_points_coast = len(coast_lon)
    num_points_isoline = len(isoline_lon)

    bounding_boxes_lonlat = []


    # crosshore wall indices
    wall_bp_file = wall_bp_file_list[bp_file_dex]
    file = open(wall_bp_file,'r')
    wall_bp_list = file.read().splitlines()
    file.close()
    wall_bp_list = [ast.literal_eval(el) for el in wall_bp_list]

    bp_file_dex += 1


    for ii in range(len(wall_bp_list)-1):
    #for ii in range(len(wall_bp_list)):

        polygon_lon = []
        polygon_lat = []
        
        polygon_lon += coast_lon[wall_bp_list[ii][0]:wall_bp_list[ii+1][0]+1]
        polygon_lat += coast_lat[wall_bp_list[ii][0]:wall_bp_list[ii+1][0]+1]
        polygon_lon += list(reversed(isoline_lon[wall_bp_list[ii][1]:wall_bp_list[ii+1][1]+1]))
        polygon_lat += list(reversed(isoline_lat[wall_bp_list[ii][1]:wall_bp_list[ii+1][1]+1]))
        polygon_lon += [coast_lon[wall_bp_list[ii][0]]]
        polygon_lat += [coast_lat[wall_bp_list[ii][0]]]
        

        bounding_boxes_lonlat.append(np.array([polygon_lon,polygon_lat]))


    print('island number: {}'.format(island_number))
    print('\n')

    #d={}
    #d["bounding_boxes_lonlat"] = bounding_boxes_lonlat
    #np.savez(bounding_boxes_file_out, **d)
    #np.save(bounding_boxes_file_out, np.array(bounding_boxes_lonlat, dtype=object), allow_pickle=True)

    file = open(bounding_boxes_file_out,'wb')
    pickle.dump(bounding_boxes_lonlat,file)
    file.close()





