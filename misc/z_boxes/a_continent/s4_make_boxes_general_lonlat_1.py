# lat/lon area calc taken from:
# https://stackoverflow.com/questions/68118907/shapely-pyproj-find-area-in-m2-of-a-polygon-created-from-latitude-and-longi


# Prevoiusly was just saving the endpoints of the "walls" running "perpendicular"
# to the coast.  But I now want to save the entire polygons of the bounding boxes,
# so that I can easily feed them into an algorithm to determine if a tracked
# particle is in one of them.

# Last version didn't store the lower wall of the polygon...

# I guess we need to "close" polygons?  i.e. copy first point as last point

# Potential bug - does this handle the case where the prevoius box's
# "closest point" is the last isoline point?  Ie there are coast points
# left, but no more isoline points?
# Wait, I think it does...

# I think that I have reversed i and j in all of this work.
# i think that i corresponds roughly to longitude, and j to lattitude

# Ok, imagine that West is up  So X = rows, Y = cols


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

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)

lon_field = d["lon_rho"]
lat_field = d["lat_rho"]
mask = d["mask_rho"]





isoline_file_in = os.path.join(output_dir,"isodistance_to_coastline_lonlat_coords_Mercator_continent.npz")

coastline_file_in = os.path.join(output_dir,"coastline_coords_Mercator_continent.npz")

bounding_boxes_file_out = os.path.join(output_dir,"bounding_boxes_lonlat_Mercator_continent.p")


# Set box area threshold (should be 300km^2, right?)
###box_area_threshold = 12.96e6
#box_area_threshold = 3.6e6 # I think this is the Raimandi size
#box_area_threshold = 10e6
box_area_threshold = 3e8
#---------------------------------------------------------------------
#---------------------------------------------------------------------

# from https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
#def PolyArea(x,y):
#    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

# For calculating area
geod = Geod(ellps="WGS84")

# Load the coastlines
d = np.load(coastline_file_in)
coastline_lonlat = d["coastline_lonlat"]

coast_lon = coastline_lonlat[:,0]
coast_lat = coastline_lonlat[:,1]

    
# Load the isolines
d = np.load(isoline_file_in)
isoline_lonlat = d["isoline_lonlat"]


isoline_lon = list(isoline_lonlat[:,0])
isoline_lat = list(isoline_lonlat[:,1])

num_points_coast = len(coast_lon)
num_points_isoline = len(isoline_lon)

bounding_boxes_lonlat = []

walk_switch = True

coast_dex = 1
isoline_dex = 1

#closees_dex = 0

while walk_switch == True:

    print("coast dex: {}, total coast points: {}".format(str(coast_dex),str(num_points_coast)))

    polygon_wall_onshore_lon = []
    polygon_wall_onshore_lat = []

    polygon_wall_onshore_lon.append(coast_lon[coast_dex-1])
    polygon_wall_onshore_lat.append(coast_lat[coast_dex-1])

    # walk up the coast
    for jj in range(coast_dex,num_points_coast):
       
        #print(coast_dex)

        # append current coast points to polygon
        polygon_wall_onshore_lon.append(coast_lon[jj])
        polygon_wall_onshore_lat.append(coast_lat[jj])

        # see which remaining (ie not already used) isoline point is nearest
        # Set a dummy "dmax" which will be replaced after the fist iteration
        dmax = 1e1000
        closest_dex = None
        for ii in range(isoline_dex,num_points_isoline):
            #dist = great_circle((coast_lat[ii],coast_lon[ii]),(isoline_lat[ii],isoline_lon[ii]))
            dist = great_circle((coast_lat[jj],coast_lon[jj]),(isoline_lat[ii],isoline_lon[ii]))
            dist = dist._Distance__kilometers
            if dist < dmax: 
                dmax = dist
                closest_dex = ii
                #print(dist)
                #print('{}  ,    {}'.format(closest_dex,dist))

        #print('{}  ,    {}'.format(closest_dex,dist))
        #walk_switch = False
        #break


        # THIS MUST CHANGE, ESPECIALLY FOR ISLANDS - WANT BALANCED BOXES, SO NEED TO 
        # CHANGE THRESHOLD AND REOCMPUTE IF LAST BOX IS TOO LARGE/SMALL

        # ALSO, CHRIS WANTS TO JUST ADD THE LAST BOX TO THE PREVIOUS, NOT MAKE A SMALL ONE
        # AS I DO HERE

        # add final polygon if we used the final isoline or coast point
        if closest_dex == num_points_isoline-1 or jj == num_points_coast-1:

            if jj == num_points_coast-1:
                polygon_lon = polygon_wall_onshore_lon[:] + isoline_lon[-1:isoline_dex-2:-1]
                polygon_lat = polygon_wall_onshore_lat[:] + isoline_lat[-1:isoline_dex-2:-1]

            else: 
                polygon_wall_onshore_lon = polygon_wall_onshore_lon + list(coast_lon[jj+1:])
                polygon_wall_onshore_lat = polygon_wall_onshore_lat + list(coast_lat[jj+1:])
                #polygon_wall_onshore_lon.append(coast_lon[jj+1:])
                #polygon_wall_onshore_lat.append(coast_lat[jj+1:])
                
                polygon_lon = polygon_wall_onshore_lon[:] + isoline_lon[-1:isoline_dex-2:-1]
                polygon_lat = polygon_wall_onshore_lat[:] + isoline_lat[-1:isoline_dex-2:-1]

            # Close polygon
            polygon_lon.append(polygon_lon[0])
            polygon_lat.append(polygon_lat[0])

            bounding_boxes_lonlat.append(np.array([polygon_lon,polygon_lat]))

            walk_switch = False

            # Print polygon area
            polygon_test_list = []
            for ii in range(len(polygon_lon)):
                polygon_test_list.append([polygon_lon[ii],polygon_lat[ii]])
            polygon_test = Polygon(polygon_test_list)
            poly_area, poly_perimeter = geod.geometry_area_perimeter(polygon_test)
            print(str(poly_area))
            
            break


        # Build a polygon to test for area
        if isoline_dex == 1:
            polygon_lon = polygon_wall_onshore_lon[:] + isoline_lon[closest_dex:isoline_dex-1:-1] + [isoline_lon[0]]
            polygon_lat = polygon_wall_onshore_lat[:] + isoline_lat[closest_dex:isoline_dex-1:-1] + [isoline_lat[0]]
        else:
            polygon_lon = polygon_wall_onshore_lon[:] + isoline_lon[closest_dex:isoline_dex-2:-1]
            polygon_lat = polygon_wall_onshore_lat[:] + isoline_lat[closest_dex:isoline_dex-2:-1]
       
        # Close polygon
        polygon_lon.append(polygon_lon[0])
        polygon_lat.append(polygon_lat[0])


        # Test polygon area
        polygon_test_list = []

        #point_list = []
    
        #line_string = LineString([Point(1, 2), Point(3, 4)])

        for ii in range(len(polygon_lon)):
            polygon_test_list.append([polygon_lon[ii],polygon_lat[ii]])
            #point_list.append(Point(polygon_lon[ii],polygon_lat[ii]))
        
        #line_string = LineString(point_list)

        polygon_test = Polygon(polygon_test_list)
        #polygon_test = Polygon(line_string)
        poly_area, poly_perimeter = geod.geometry_area_perimeter(polygon_test)
        
        #print(poly_area)

        if poly_area >= box_area_threshold:

            bounding_boxes_lonlat.append(np.array([polygon_lon,polygon_lat]))

            isoline_dex = closest_dex + 1
            coast_dex = jj + 1

            #print(str(poly_area))
            print('polygon area: {},  threshold: {}'.format(poly_area,box_area_threshold))
            print('\n')

            break


file = open(bounding_boxes_file_out,'wb')
#pickle.dump(walls_ij,file)
pickle.dump(bounding_boxes_lonlat,file)
file.close()






