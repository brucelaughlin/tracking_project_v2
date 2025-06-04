
from netCDF4 import Dataset

rg = Dataset("test.nc", "a", format="NETCDF3_CLASSIC")

lat = rg.createDimension("lat", 100)
lon = rg.createDimension("lon", 100)



#rg.close()
