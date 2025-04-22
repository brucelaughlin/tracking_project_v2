import netCDF4
import numpy as np
import matplotlib.pyplot as plt


#land_type = 'continent'
#land_type = 'islands'

output_dir = 'z_output/'

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)
lon_rho = d["lon_rho"]
lat_rho = d["lat_rho"]
mask_rho = d["mask_rho"]


coastline_file_in = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_continent/z_output/coastline_coords_psi_file_wc15n_continent.npz"
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
ax.pcolormesh(lon_rho,lat_rho,mask_rho,shading="nearest")

ax.plot(coastline_lonlat_continent[:,0],coastline_lonlat_continent[:,1])
ax.plot(coastline_lonlat_island_1[:,0],coastline_lonlat_island_1[:,1])
ax.plot(coastline_lonlat_island_2[:,0],coastline_lonlat_island_2[:,1])
ax.plot(coastline_lonlat_island_3[:,0],coastline_lonlat_island_3[:,1])

ax.axis('image')
plt.show()









