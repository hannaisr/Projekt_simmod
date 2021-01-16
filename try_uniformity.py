import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x=[]
y=[]
z=[]

alpha = 0.1

# Old version
# for i in range(1000):
#     theta = rnd.uniform(0,np.pi)
#     phi = rnd.uniform(0,2*np.pi)
#     x.append(np.sin(theta)*np.cos(phi))
#     y.append(np.sin(theta)*np.sin(phi))
#     z.append(np.cos(theta))

# New try
for i in range(1000):
    theta = 2*np.pi*rnd.uniform(0,1)
    phi = np.arccos(2*rnd.uniform(0,1)-1)
    x.append(np.sin(theta)*np.cos(phi))
    y.append(np.sin(theta)*np.sin(phi))
    z.append(np.cos(theta))

# Uniform distribuation?
# for i in range(1000):
#     x_ny = rnd.uniform(-1,1)
#     y_ny = rnd.uniform(-1,1)
#     z_ny = rnd.uniform(np.cos(alpha),1)
#     norm = np.sqrt(x_ny**2+y_ny**2+z_ny**2)
#     x.append(x_ny/norm)
#     y.append(y_ny/norm)
#     z.append(z_ny/norm)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z,'.')
# ax.plot(x_back,y_back,z_back,'*')
ax.plot(0,0,0,'*') # origo
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')
plt.show()
