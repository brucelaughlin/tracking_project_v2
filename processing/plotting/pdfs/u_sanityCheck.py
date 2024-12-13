import matplotlib.pyplot as plt
import numpy as np


z = np.ones((3,6))
x = np.arange(-.5,6,1)
y = np.arange(-.5,3,1)


z[1,:] *= 2
z[2,:] *= 3

#fig,ax = plt.subplots()
###blah = ax.pcolormesh(x,y,z)
#bah = ax.pcolormesh(y,x,z.T)
#fig.colorbar(blah,ax=ax)
#plt.show()

row_sums = z.sum(axis=1)
col_sums = z.sum(axis=0)

print(row_sums)
print("\n")
print(col_sums)


