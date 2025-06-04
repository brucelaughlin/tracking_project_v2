import netCDF4
import numpy as np
import matplotlib.pyplot as plt


#land_type = 'continent'
#land_type = 'islands'

output_dir = 'z_output/'

#grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"
grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_landNAN.npz"
#grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid_noInland.npz"

d = np.load(grid_file)
lon_rho = d["lon_rho"]
lat_rho = d["lat_rho"]
mask_rho = d["mask_rho"]


a = plt.pcolormesh(mask_rho)
#pcolormesh(lon_rho,lat_rho,mask_rho,shading="nearest")

plt.colorbar(a)

plt.show()









