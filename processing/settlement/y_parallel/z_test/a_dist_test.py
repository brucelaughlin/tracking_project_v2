
import numpy as np
from geopy.distance import geodesic

t_1 = np.array([
    [37.7749, -122.4194],  
    [34.0522, -118.2437],  
    [32.7157, -96.6203],  
])

t_2 = np.copy(t_1)

t_2[1] += 5
t_2[2] += 10


total_distance = np.zeros(len(t_1))

for ii in range(len(t_1)):
    #for jj in range(np.shape(t_1)[1]):
    point1 = t_1[ii]
    point2 = t_2[ii]
    distance = geodesic(point1, point2).km  # Returns distance in kilometers
    total_distance[ii] += distance


#print(f"Total trajectory distance: {total_distance:.2f} km")
