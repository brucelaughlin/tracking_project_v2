
import os
from pathlib import Path
import pickle
import numpy as np
import matplotlib.pyplot as plt


#-------------------- EDIT THESE -------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

cwd = os.getcwd()
output_dir = os.path.join(str(cwd),"z_output")
Path(output_dir).mkdir(parents=True, exist_ok=True)

grid_file = "/home/blaughli/tracking_project_v2/misc/z_boxes/z_output/mercator_diy_grid.npz"

d = np.load(grid_file)

lon_field = d["lon_rho"]
lat_field = d["lat_rho"]
mask = d["mask_rho"]



bounding_boxes_file_in = "/home/blaughli/tracking_project_v2/misc/z_boxes/a_continent/z_output/bounding_boxes_lonlat_Mercator_continent.p" 

#---------------------------------------------------------------------
#---------------------------------------------------------------------


# Load the boxes
file = open(bounding_boxes_file_in,'rb')
boxes_lonlat = pickle.load(file)
file.close

fig, ax = plt.subplots()
ax.pcolormesh(lon_field,lat_field,mask,shading="nearest")

for box in boxes_lonlat:
    if box is not None:
       #ax.plot(box[1],box[0])
       ax.plot(box[0],box[1])

ax.axis('image')
plt.show()









