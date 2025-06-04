# Data is regular!  So can use regular grid interpolator


from scipy.interpolate import RegularGridInterpolator
#from scipy.interpolate import LinearNDInterpolator
import netCDF4
import numpy as np
import os
import matplotlib.pyplot as plt

input_model_nc_file = "/home/blaughli/symbolic_links_ROMS/z_testing/Mercator/global-reanalysis-phy-001-030-daily_1993_addAttr.nc"
#input_model_nc_file = "/home/blaughli/symbolic_links_ROMS/z_testing/Mercator/global-reanalysis-phy-001-030-daily_1993.nc"

input_bathymetry_nc_file = "/home/blaughli/tracking_project_v2/misc/z_modifiers/z_testing/gebco_2024_n61.2714_s14.9579_w-151.3085_e-104.8933.nc"

output_model_depths_file = "/home/blaughli/tracking_project_v2/misc/z_modifiers/z_testing/model_seafloor_depths_from_gebco_2024.npz"

#with netCDF4.Dataset(input_bathymetry_nc_file,'r',format="NETCDF3_CLASSIC") as ncfile:
#    bathymetry = ncfile["elevation"]

#ncfile = netCDF4.Dataset(input_model_nc_file,'r',format="NETCDF3_CLASSIC")
with netCDF4.Dataset(input_model_nc_file,'r',format="NETCDF3_CLASSIC") as ncfile:
    land_binary_mask = np.array(ncfile["land_binary_mask"])
    lon_model_1D = np.array(ncfile["longitude"])
    lat_model_1D = np.array(ncfile["latitude"])
#ncfile.close()

lon_model,lat_model = np.meshgrid(lon_model_1D,lat_model_1D)

#a=plt.pcolormesh(lon_model,lat_model,land_binary_mask)
#plt.colorbar(a)
#plt.show()


#ncfile = netCDF4.Dataset(input_bathymetry_nc_file,'r',format="NETCDF3_CLASSIC")
with netCDF4.Dataset(input_bathymetry_nc_file,'r',format="NETCDF3_CLASSIC") as ncfile:
    bathymetry = np.array(ncfile["elevation"])
    lon_bathy_1D = np.array(ncfile["lon"])
    lat_bathy_1D = np.array(ncfile["lat"])
#ncfile.close()

lon_bathy,lat_bathy = np.meshgrid(lon_bathy_1D,lat_bathy_1D)

bathymetry[bathymetry > 0] = 0
#bathymetry[bathymetry > 0] = np.nan
bathymetry *= -1


#fig,ax = plt.subplots(2)
#a = ax[0].pcolormesh(lon_bathy,lat_bathy,bathymetry)
#plt.colorbar(a)
#b = ax[1].pcolormesh(lon_model,lat_model,land_binary_mask)
#plt.colorbar(b)
#plt.show()

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"
d = np.load(grid_file)
lon_grid = d["lon_rho"]
lat_grid = d["lat_rho"]
mask_grid = d["mask_rho"]

#b = ax[1].pcolormesh(lon_grid,lat_grid,mask_grid)
#plt.show()


interp = RegularGridInterpolator((lat_bathy_1D,lon_bathy_1D), bathymetry)
#interp = RegularGridInterpolator((lon_bathy_1D,lat_bathy_1D), bathymetry)
#interp = RegularGridInterpolator((lon_bathy.flatten(),lat_bathy.flatten()), bathymetry)
#interp = RegularGridInterpolator(list(zip(lon_bathy.flatten(),lat_bathy.flatten())), bathymetry)

#model_grid_depths = interp(lat_bathy,lon_bathy)
model_grid_depths = interp(np.array((lat_model.flatten(),lon_model.flatten())).T)

model_grid_depths=np.reshape(model_grid_depths,(np.shape(lon_model)[0],np.shape(lon_model)[1]))


d = {}
d["model_grid_depths"] = model_grid_depths
np.savez(output_model_depths_file, **d)

# THE GOOD FINAL PLOT!
#fig,ax = plt.subplots(2)
#a = ax[0].pcolormesh(lon_bathy,lat_bathy,bathymetry)
#ax[0].set_title("external bathymetry data")
#plt.colorbar(a)
#b = ax[1].pcolormesh(lon_model,lat_model,model_grid_depths)
#ax[1].set_title("bathymetry interpolated to model grid")
#plt.colorbar(b)
#plt.show()





###interp = LinearNDInterpolator(list(zip(lon_bathy.flatten(),lat_bathy.flatten())), bathymetry)
#####interp = LinearNDInterpolator(list(zip(lon_bathy.flatten(),lat_bathy.flatten())), bathymetry.flatten())




#with netCDF4.Dataset(input_model_nc_file,'a',format="NETCDF3_CLASSIC") as ncfile:
#    land_binary_mask = ncfile.createVariable("land_binary_mask","b",("latitude","longitude",))
#    land_binary_mask[:] = mask
