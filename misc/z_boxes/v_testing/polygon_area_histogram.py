
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

geod = Geod(ellps="WGS84")

polygon_csv_file_path = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/Mercator_300km2_settling_polygons.txt"

square_m_per_square_km = 1e6

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

list_of_polygon_areas = []

first_island_index = 0

current_polygon_mean_lat = np.mean(list_of_polygon_vertex_lonlat_arrays[0][:,1])

for jj in range(len(list_of_polygon_vertex_lonlat_arrays)):
    
    current_polygon = list_of_polygon_vertex_lonlat_arrays[jj]

    if np.mean(list_of_polygon_vertex_lonlat_arrays[jj][:,1]) < current_polygon_mean_lat - 5:
        first_island_index = jj
    
    current_polygon_mean_lat = np.mean(list_of_polygon_vertex_lonlat_arrays[jj][:,1])

    polygon_test = Polygon(current_polygon)
    poly_area, poly_perimeter = geod.geometry_area_perimeter(polygon_test)
    list_of_polygon_areas.append(poly_area/square_m_per_square_km)

    if poly_area > 500 or poly_area < 250:
        print(jj)
        print(poly_area)


# Filter out island polygons, since those are currently made by hand (as opposed to the algorithm)
list_of_polygon_areas_continent = list_of_polygon_areas[0:first_island_index + 1]

#np.histogram(list_of_polygon_areas)
hist, bin_edges = np.histogram(list_of_polygon_areas_continent)
###hist, bin_edges = np.histogram(list_of_polygon_areas)

csv_file_path_stem = Path(polygon_csv_file_path).stem

plt.hist(list_of_polygon_areas_continent, bins=40)
###plt.hist(list_of_polygon_areas_continent, bins='auto')
plt.title(f"Histogram of continental polygon areas\nfile: {csv_file_path_stem}")
###plt.title(f"Histogram of polygon areas\nfile: {csv_file_path_stem}")

plt.xlabel("Area (square kilometers)")
plt.ylabel(f"Count ({len(list_of_polygon_areas)} total polygons)")

plt.show()



