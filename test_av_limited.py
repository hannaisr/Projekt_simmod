import numpy as np
import random as rnd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# alpha = 0.1*np.pi
theta0 = rnd.uniform(0,np.pi)
phi0 = rnd.uniform(0,2*np.pi)
# x0 = np.sin(theta0)*np.cos(phi0)
# y0 = np.sin(theta0)*np.sin(phi0)
# z0 = np.cos(theta0)

rho = 0.5
r = 1

theta_back = np.pi-theta0
phi_back = np.pi+phi0

x_back = np.sin(theta_back)*np.cos(phi_back)
y_back = np.sin(theta_back)*np.sin(phi_back)
z_back = np.cos(theta_back)

#In local coordinate system with "going back" vector as z axis:
#Generate new bead center
alpha = np.arccos((r**2 + r**2 - (2*rho) ** 2) / (2*r*r))
theta = rnd.uniform(alpha, np.pi)
phi = rnd.uniform(0, 2*np.pi)

#Translate bead center position into global coordinate system.
# First: define the local coordinate system in terms of the global
z_hat = [np.sin(theta_back) * np.cos(phi_back), np.sin(theta_back) * np.sin(phi_back),np.cos(theta_back)]  # vector between the most recent two beads
y_hat = [-np.sin(phi_back), np.cos(phi_back), 0] # changed x->y
# x_hat = [x_hat[i]* 1 / np.sin(theta_back) for i in range(3)]  # Orthogonal to z_hat
x_hat = [-np.cos(theta_back) * np.cos(phi_back),-np.cos(theta_back) * np.sin(phi_back), np.sin(theta_back)] # Orthogonal to z_hat, x_hat #changed y->x
#    y_hat = [y_hat[i]* 1 / np.sin(theta_back) for i in range(3)]  # Orthogonal to z_hat

theta_list = []
phi_list=[]
for i in range(5000):
	theta_list.append(rnd.uniform(alpha,np.pi))
	phi_list.append(rnd.uniform(0,2*np.pi))

# print(theta_list,phi_list)
#print(sum([y_hat[i] * y_hat[i] for i in range(3)]))
#Second: project the bead center position onto the global axes and translate it to origo
x=[]
y=[]
z=[]
for i in range(len(theta_list)):
	x.append(r*(np.cos(theta_list[i])*z_hat[0]+np.sin(theta_list[i])*np.cos(phi_list[i])*y_hat[0]+np.sin(theta_list[i])*np.sin(phi_list[i])*x_hat[0]))
	y.append(r*(np.cos(theta_list[i])*z_hat[1]+np.sin(theta_list[i])*np.cos(phi_list[i])*y_hat[1]+np.sin(theta_list[i])*np.sin(phi_list[i])*x_hat[1]))
	z.append(r*(np.cos(theta_list[i])*z_hat[2]+np.sin(theta_list[i])*np.cos(phi_list[i])*y_hat[2]+np.sin(theta_list[i])*np.sin(phi_list[i])*x_hat[2]))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z,'.')
ax.plot(x_back,y_back,z_back,'*')
ax.plot(0,0,0,'*') # origo
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')

plt.show()
