import netCDF4
import numpy as np
import matplotlib.pyplot as plt

output_dir = 'z_output/'

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

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
m =ax.pcolormesh(lon_rho,lat_rho,mask_rho_new,shading="nearest")

ax.plot(coastline_lonlat_continent[:,0],coastline_lonlat_continent[:,1],c='k')
ax.plot(coastline_lonlat_island_1[:,0],coastline_lonlat_island_1[:,1],c='k')
ax.plot(coastline_lonlat_island_2[:,0],coastline_lonlat_island_2[:,1],c='k')
ax.plot(coastline_lonlat_island_3[:,0],coastline_lonlat_island_3[:,1],c='k')

pt = 100
ax.scatter(lon_rho[pt,pt],lat_rho[pt,pt],c='k')

ax.scatter(lon_psi[pt,pt],lat_psi[pt,pt],c='r')
ax.scatter(lon_psi[pt,pt-1],lat_psi[pt,pt-1],c='c')
ax.scatter(lon_psi[pt-1,pt-1],lat_psi[pt-1,pt-1],c='g')
ax.scatter(lon_psi[pt-1,pt],lat_psi[pt-1,pt],c='b')


ax.axis('image')

plt.colorbar(m)

plt.show()


