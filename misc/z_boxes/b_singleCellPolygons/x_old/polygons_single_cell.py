import netCDF4
import numpy as np
import matplotlib.pyplot as plt

output_dir = 'z_output/'

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_noInland.npz"
grid_file_plot = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)
lon_rho = d["lon_rho"]
lat_rho = d["lat_rho"]
mask_rho = d["mask_rho"]
mask_psi = d["mask_psi"]
lon_psi = d["lon_psi"]
lat_psi = d["lat_psi"]

d = np.load(grid_file_plot)
mask_rho_plot = d["mask_rho"]

mask_rho_plot[mask_rho_plot == 0] = np.nan

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
#mask_rho_new[mask_rho == 0] = 15
#mask_rho_new[mask_rho_new == 10] = 15


#first attempt at making single-cell polygons; this doesn't do ordering in any intelligent way

list_of_polygon_vertex_lonlat_lists = []

ii_max = 394
ii_min = 145

for ii in range(ii_min,ii_max):
#for ii in range(np.shape(mask_rho_new)[0]):
    for jj in range(1,np.shape(mask_rho_new)[1]-1):
    #for jj in range(np.shape(mask_rho_new)[1]):

        if np.logical_not(np.isnan(mask_rho_new[ii,jj])):

            polygon_lon = []
            polygon_lat = []
            
            polygon_lon.append(lon_psi[ii,jj])
            polygon_lon.append(lon_psi[ii,jj-1])
            polygon_lon.append(lon_psi[ii-1,jj-1])
            polygon_lon.append(lon_psi[ii-1,jj])
            polygon_lon.append(lon_psi[ii,jj])
            
            polygon_lat.append(lat_psi[ii,jj])
            polygon_lat.append(lat_psi[ii,jj-1])
            polygon_lat.append(lat_psi[ii-1,jj-1])
            polygon_lat.append(lat_psi[ii-1,jj])
            polygon_lat.append(lat_psi[ii,jj])
            
            #list_of_polygon_vertex_longitudes.append(lon_psi[ii,jj],lat_psi[ii,jj],c='r')
            #list_of_polygon_vertex_longitudes.append(lon_psi[ii,jj-1],lat_psi[ii,jj-1],c='c')
            #list_of_polygon_vertex_longitudes.append(lon_psi[ii-1,jj-1],lat_psi[ii-1,jj-1],c='g')
            #list_of_polygon_vertex_longitudes.append(lon_psi[ii-1,jj],lat_psi[ii-1,jj],c='b')

            list_of_polygon_vertex_lonlat_lists.append(np.array([polygon_lon,polygon_lat]))

#            print(ii)
#            print(jj)
#            break
#
#    else:
#        continue
#    break


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
m = ax.pcolormesh(lon_rho,lat_rho,mask_rho_plot,shading="nearest")
#m = ax.pcolormesh(lon_rho,lat_rho,mask_rho,shading="nearest")
#m = ax.pcolormesh(lon_rho,lat_rho,mask_rho_new,shading="nearest")

###ax.plot(coastline_lonlat_continent[:,0],coastline_lonlat_continent[:,1],c='k')
###ax.plot(coastline_lonlat_island_1[:,0],coastline_lonlat_island_1[:,1],c='k')
###ax.plot(coastline_lonlat_island_2[:,0],coastline_lonlat_island_2[:,1],c='k')
###ax.plot(coastline_lonlat_island_3[:,0],coastline_lonlat_island_3[:,1],c='k')

for ii in range(len(list_of_polygon_vertex_lonlat_lists)):
    ax.plot(list_of_polygon_vertex_lonlat_lists[ii][0,:],list_of_polygon_vertex_lonlat_lists[ii][1,:], c='y')
    #ax.plot(list_of_polygon_vertex_lonlat_lists[ii][0,:],list_of_polygon_vertex_lonlat_lists[ii][1,:], c='k')

#ax.plot(list_of_polygon_vertex_lonlat_lists[0][0,:],list_of_polygon_vertex_lonlat_lists[0][1,:], c='k')

ax.axis('image')

#plt.colorbar(m)

plt.show()


