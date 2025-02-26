import numpy as np
import matplotlib.pyplot as plt

z_array = np.arange(0,100)
#z_array = np.arange(0,100)
target_depth = 30
#vkick_max = 0.05
#random_vkick_max = 0.005
timestep_seconds = 60*60

vkick_max = 0.005
random_vkick_max = vkick_max/10


# In the model I multiply by -1... maybe bad!
#z_array *= -1

tanh_velocities = np.tanh(np.pi/target_depth * (z_array - target_depth)) * vkick_max

line_velocities = (z_array - target_depth)/timestep_seconds
#line_velocities = -1 * (z_array - target_depth)/timestep_seconds


vel_use = np.minimum(abs(tanh_velocities), abs(line_velocities))

#idx = np.argwhere(np.diff(np.sign(tanh_velocities - line_velocities))).flatten()
#plt.plot(z_array[idx], line_velocities[idx], 'ro')
#plt.plot(z_array[idx], tanh_velocities[idx], 'ro')


plt.plot(tanh_velocities)
plt.plot(line_velocities)

sign_mask = np.ones(len(z_array))
sign_mask[z_array < target_depth] *= -1

#plt.plot(abs(tanh_velocities))
#plt.plot(abs(line_velocities))

plt.plot(vel_use * sign_mask)

plt.show()



