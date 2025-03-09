import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("vswimmax", type=float)
parser.add_argument("targetdepth", type=float)
args = parser.parse_args()

vswim_max = args.vswimmax
target_depth = args.targetdepth

swim_cut_depth_distance_from_target_depth = 5

upper_swim_cut_depth = target_depth - swim_cut_depth_distance_from_target_depth
lower_swim_cut_depth = target_depth + swim_cut_depth_distance_from_target_depth

noise_region_radius = 5
noise_std = noise_region_radius/4
#noise_std = noise_region_radius/2


z_array = np.arange(0,100)
timestep_seconds = 60*60


sign_mask = np.ones(len(z_array))
sign_mask[z_array < target_depth] *= -1


swim_tanh = vswim_max * (-1 * (1 - (1 + np.tanh(z_array - upper_swim_cut_depth))/2  - (1 + np.tanh(z_array - lower_swim_cut_depth))/2))
swim_line = (z_array - target_depth)/timestep_seconds

vel_use = np.minimum(abs(swim_tanh), abs(swim_line))

swim_noise = target_depth - np.random.normal(target_depth, noise_std/timestep_seconds, len(z_array))

vel_total = vel_use * sign_mask + swim_noise


plt.plot(swim_tanh)
plt.plot(swim_line)

#plt.plot(vel_total)
plt.plot(vel_use * sign_mask)

plt.show()



