# Need to create fields for "dist2coast.m", which I previously got from netcdf grid files

import scipy
import numpy as np
import argparse
from pathlib import Path
import os
import netCDF4

parser = argparse.ArgumentParser()
parser.add_argument("gridfile", type=str)
args = parser.parse_args()

grid_file = args.gridfile


output_dir = "z_output/"
Path(output_dir).mkdir(parents=True, exist_ok=True)

output_grid_file = os.path.join(output_dir,f"{Path(grid_file).stem}.mat")


if os.path.splitext(grid_file)[-1] == ".npz":
    d = np.load(grid_file)
    lon_rho = d["lon_rho"]
    lat_rho = d["lat_rho"]
    mask_rho = d["mask_rho"]
elif os.path.splitext(grid_file)[-1] == ".nc":
    d = netCDF4.Dataset(grid_file)
    lon_rho = np.array(d["lon_rho"])
    lat_rho = np.array(d["lat_rho"])
    mask_rho = np.array(d["mask_rho"])
    d.close()
else:
    print("Grid file is expected to have either '.nc' or '.npz' as its extension")
    sys.exit(1)


d = {}
d["mask_rho"] = mask_rho
d["lon_rho"] = lon_rho
d["lat_rho"] = lat_rho

scipy.io.savemat(output_grid_file, d)
