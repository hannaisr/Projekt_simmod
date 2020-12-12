import random as rnd
import numpy as np
import matplotlib.pyplot as plt # For 3D plot
from mpl_toolkits.mplot3d import Axes3D # For 3D plot

class Walker():
    """Walked positions are stored in three lists, one for each coordinate."""
    # Initial position
    x0 = 0
    y0 = 0
    z0 = 0

    # Define step length
    r = 1

    # Initiate the list storing the coordinates of the walk
    x, y, z = [x0], [y0], [z0]

    def walk_one_step(self):
        """Implemented by child class"""
        pass

    def non_avoiding_walk(self,nsteps=100):
        """Walk nsteps steps of non-self-avoiding random walk."""
        for i in range(nsteps):
            self.walk_one_step()

    def self_avoiding_walk(self,nsteps=100):
        """Implemented by child class"""
        pass

    def restart(self):
        self.x, self.y, self.z = [self.x0], [self.y0], [self.z0]

    def plot_the_walk(self):
        """Plots the walk in a 3D line plot."""
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.x,self.y,self.z)
        plt.show()

class Grid_walker(Walker):
    def walk_one_step(self):
        # Append the same coordinates as last step to coordinate list
        self.x.append(self.x[-1]), self.y.append(self.y[-1]), self.z.append(self.z[-1])
        # Get walking direction
        direction = rnd.randint(0,5)
        # Update the coordinates
        if direction == 0:
            self.x[-1] += 1
        elif direction == 1:
            self.x[-1] -= 1
        elif direction == 2:
            self.y[-1] += 1
        elif direction == 3:
            self.y[-1] -= 1
        elif direction == 4:
            self.z[-1] += 1
        elif direction == 5:
            self.z[-1] -= 1

    def self_avoiding_walk(self,nsteps=100):
        """Walk nsteps steps of self-avoiding random walk"""
        pass

class Freely_jointed_chain(Walker):
    # Define radius of spheres at ends for self-avoiding walk (could perhaps be put in __init__())
    R = 1/2

    def walk_one_step(self):
        # Append the same coordinates as last step to coordinate list
        self.x.append(self.x[-1]), self.y.append(self.y[-1]), self.z.append(self.z[-1])
        # Get walking direction
        theta = rnd.uniform(0,np.pi)
        phi = rnd.uniform(0,2*np.pi)
        # Update the coordinates
        self.x[-1] += self.r*np.sin(theta)*np.cos(phi)
        self.y[-1] += self.r*np.sin(theta)*np.sin(phi)
        self.z[-1] += self.r*np.cos(theta)

    def self_avoiding_walk(self,nsteps=100):
        """Walk nsteps steps of self-avoiding random walk"""
        pass

# gridwalk = Grid_walker()
# gridwalk.non_avoiding_walk(nsteps=500)
# gridwalk.plot_the_walk()

chainwalk = Freely_jointed_chain()
chainwalk.non_avoiding_walk(nsteps=500)
chainwalk.plot_the_walk()
