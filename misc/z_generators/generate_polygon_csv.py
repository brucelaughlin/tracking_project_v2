# Write polygons lat/lon to a csv file

import numpy as np
import os
import pickle
import csv

output_file = 'wc15n_300km2_settling_polygons.txt'

input_file_dir = '/home/blaughli/tracking_project_v2/input_files'

input_file_list = []
input_file_list.append('bounding_boxes_lonlat_coords_psi_coastline_wc15n_continent.p')
input_file_list.append('bounding_boxes_lonlat_wc15n_island_number_4.p')
input_file_list.append('bounding_boxes_lonlat_wc15n_island_number_5.p')
input_file_list.append('bounding_boxes_lonlat_wc15n_island_number_6.p')
input_file_list.append('bounding_boxes_lonlat_wc15n_island_number_7.p')
input_file_list.append('bounding_boxes_lonlat_wc15n_island_number_8.p')


field_names = ['cell #', ' vertex #', ' lat', ' lon']

csv_data = []

csv_data.append(field_names)

cell_number = 0



for input_file in input_file_list:

    polygons_file_in = os.path.join(input_file_dir, input_file)

    file = open(polygons_file_in,'rb')
    polygons_lonlat = pickle.load(file)
    file.close

    for polygon_lonlat in polygons_lonlat:
        
        cell_number += 1

        for vertex_number in range(np.shape(polygon_lonlat)[1]):   
        #for vertex_number in range(np.shape(polygon_lonlat)[1] - 1):   

            vertex_number_print = vertex_number + 1

            data = []
            data.append(f"{cell_number:03d}")
            data.append(f" {vertex_number_print:02d}")
            data.append(f" {polygon_lonlat[1,vertex_number]:.6f}")
            data.append(f" {polygon_lonlat[0,vertex_number]:.6f}")

            csv_data.append(data)


with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(csv_data)




