# just check that columns do actually add to 1

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as tkr


#a=np.random.rand(3,3)
a=np.array([[5,25,70],[89,1,10],[10,1,89]])
#a=np.array([[5,25,70],[7.5,17.5,75],[10,1,89]])

row_sums = a.sum(axis=1)
a = a / row_sums
#a = a / row_sums[:, np.newaxis]

#a=np.array([[.05,.25,.7],[.075,.175,.75],[.1,.01,.89]])
##a=np.array([[.1,.3,.6],[.1,.3,.6],[.1,.3,.6]])

print("test data:")
print(a)

print("test data row sums:")
print(np.sum(a,axis=1))

plt.pcolormesh(a.T)
plt.colorbar()
plt.show()


