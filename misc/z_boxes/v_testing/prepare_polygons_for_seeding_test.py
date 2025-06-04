
import pdb
import time

import datetime
import netCDF4
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as plt_path
from scipy.interpolate import interp1d
from geopy.distance import great_circle
import scipy.interpolate as spint
from os import listdir
from os.path import isfile, join
import sys
import argparse
from pathlib import Path
import os
import math

script_time_start = time.time()

#---------------------------------------------------------------------
#---------------------------------------------------------------------
# PARAMETERS and CONSTANTS

# Opendrift swaps the particle and time dimensions relative to the tracking code that Patrick used


parser = argparse.ArgumentParser()
parser.add_argument("--polygoncsvfilepath",type=str, help='blah')
args = parser.parse_args()

polygon_csv_file_path = args.polygoncsvfilepath


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
num_polygons = len(list_of_polygon_vertex_lonlat_arrays)
#---------------------------------------------------------------------



