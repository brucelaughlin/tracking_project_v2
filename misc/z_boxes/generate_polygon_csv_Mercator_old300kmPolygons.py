# Write polygons lat/lon to a csv file

import numpy as np
import os
import pickle
import csv

output_dir = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output"

output_file = os.path.join(output_dir,'Mercator_300km2_settling_polygons.txt')


input_file_list = []
input_file_list.append('/home/blaughli/tracking_project_v2/misc/z_boxes/a_continent/z_output/bounding_boxes_lonlat_Mercator_continent.p')
input_file_list.append('/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/z_output/bounding_boxes_lonlat_Mercator_island_number_2.p')
input_file_list.append('/home/blaughli/tracking_project_v2/misc/z_boxes/a_islands/z_output/bounding_boxes_lonlat_Mercator_island_number_3.p')

field_names = ['cell #', ' vertex #', ' lat', ' lon']

csv_data = []

csv_data.append(field_names)

cell_number = 0



for polygons_file_in in input_file_list:
#for input_file in input_file_list:

    #polygons_file_in = os.path.join(input_file_dir, input_file)

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




