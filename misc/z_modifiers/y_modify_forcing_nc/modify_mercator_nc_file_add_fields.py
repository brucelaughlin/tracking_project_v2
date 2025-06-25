
import netCDF4
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("inputmodelforcingfile", type=str)

args = parser.parse_args()
input_model_forcing_file = args.inputmodelforcingfile
#input_model_forcing_file = "/home/blaughli/symbolic_links_ROMS/z_testing/Mercator/global-reanalysis-phy-001-030-daily_1993_addAttr.nc"


#input_model_grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"
input_model_depths_file = "/home/blaughli/tracking_project_v2/misc/z_modifiers/y_modify_forcing_nc/z_output/model_seafloor_depths_Mercator.npz"


# --------------------------------------------------------------
#d = np.load(input_model_grid_file)
#mask_rho = d["mask_rho"]
# --------------------------------------------------------------
# Opendrift seems to want ocean points to have a mask value of 0, and land a value of 1.  So, for now, I need to reverse those in my custom mask
#dex_ocean = mask_rho == 1
#dex_land = mask_rho == 0
#mask_rho[dex_ocean] = 0
#mask_rho[dex_land] = 1
# --------------------------------------------------------------

d = np.load(input_model_depths_file)
model_grid_depths = d["model_grid_depths"] 
mask_rho = d["mask_rho"]

with netCDF4.Dataset(input_model_forcing_file,'a',format="NETCDF3_CLASSIC") as ncfile:
    land_binary_mask = ncfile.createVariable("land_binary_mask","b",("latitude","longitude",))
    land_binary_mask[:] = mask_rho
    land_binary_mask.standard_name = "land_binary_mask"

    sea_floor_depth_below_sea_level = ncfile.createVariable("sea_floor_depth_below_sea_level","f",("latitude","longitude",))
    sea_floor_depth_below_sea_level[:] = model_grid_depths
    sea_floor_depth_below_sea_level.standard_name = "sea_floor_depth_below_sea_level"
