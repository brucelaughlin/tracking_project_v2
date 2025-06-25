
import matplotlib.pyplot as plt
import numpy as np
import argparse
import netCDF4

parser = argparse.ArgumentParser()
parser.add_argument("gridfile", type=str)
parser.add_argument("latlonswitch", type=int)
parser.add_argument("nonromsswitch", nargs="?")
args = parser.parse_args()

grid_file = args.gridfile
latlon_switch = bool(args.latlonswitch)
non_roms_switch = args.nonromsswitch

if non_roms_switch is not None:
    d = np.load(grid_file)
    lat = d["lat_rho"]
    lon = d["lon_rho"]
    mask = d["mask_rho"]
else:
    d = netCDF4.Dataset(grid_file,"r")
    lat = np.array(d["lat_rho"])
    lon = np.array(d["lon_rho"])
    mask = np.array(d["mask_rho"])
    d.close()

if latlon_switch:
    plt.pcolormesh(lon,lat,mask)
else:
    plt.pcolormesh(mask)

plt.show()

