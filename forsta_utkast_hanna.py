import random as rnd
import numpy as np
import matplotlib.pyplot as plt # For 3D plot
from mpl_toolkits.mplot3d import Axes3D # For 3D plot

class Walker():
    """Walked positions are stored in a list"""
    # Initiate list of visited points
    origin = (0.,0.,0.)
    visited_points = [origin]

    # Define step length
    r = 1.

    def walk_one_step(self):
        """Implemented by child class"""
        pass

    def walk_without_avoid(self,nsteps=100):
        """Walk nsteps steps of non-self-avoiding random walk."""
        for i in range(nsteps):
            self.walk_one_step()

    def walk_with_self_avoid(self,nsteps=100):
        """Implemented by child class"""
        pass

    def restart(self):
        """Resets list of visited points."""
        self.visited_points = [self.origin]

    def plot_the_walk(self):
        """Plots the walk in 3D."""
        visited_points = np.array(self.visited_points)
        x = [i[0] for i in visited_points]
        y = [i[1] for i in visited_points]
        z = [i[2] for i in visited_points]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x,y,z)
        plt.show()

    def get_end_to_end_distance(self):
        """Calculate end-to-end distance of already walked walk."""
        last_point = np.array(self.visited_points[-1][:])
        return np.sqrt((last_point[0]-self.origin[0])**2+(last_point[1]-self.origin[1])**2+(last_point[2]-self.origin[2])**2)

    def get_multiple_end_to_end_distances(self,nwalks=10):
        """Returns a list of end-to-end distances for nwalks number of walks"""
        # TODO implement for self avoiding walk as well
        etedist_list = np.zeros(nwalks)
        for i in range(nwalks):
            self.restart()
            self.walk_without_avoid()
            etedist_list[i] = self.get_end_to_end_distance()
        return etedist_list

class Grid_walker(Walker):
    def walk_one_step(self):
        current_pos = np.array(self.visited_points[-1])
        # Get walking direction
        direction = rnd.randint(0,5)
        # Update the coordinates
        if direction == 0:
            current_pos[0] += 1
        elif direction == 1:
            current_pos[0] -= 1
        elif direction == 2:
            current_pos[1] += 1
        elif direction == 3:
            current_pos[1] -= 1
        elif direction == 4:
            current_pos[2] += 1
        elif direction == 5:
            current_pos[2] -= 1
        # Update list of visited points
        self.visited_points.append(tuple(current_pos))

    def walk_with_self_avoid(self,nsteps=100):
        """Walk nsteps steps of self-avoiding random walk"""
        for i in range(nsteps):
            self.walk_one_step()
            current_pos = self.visited_points[-1]
            if any(t == current_pos for t in self.visited_points[:-1]):
                print('Managed to walk',len(self.visited_points)-1,'steps')
                break

class Freely_jointed_chain(Walker):
    # Define radius of spheres at ends for self-avoiding walk (could perhaps be put in __init__())
    R = 1/2

    def walk_one_step(self):
        current_pos = np.array(self.visited_points[-1])
        # Get walking direction
        theta = rnd.uniform(0,np.pi)
        phi = rnd.uniform(0,2*np.pi)
        print('theta:',theta,'phi:',phi,'r:',self.r)
        # Update the coordinates
        current_pos[0] += self.r*np.sin(theta)*np.cos(phi)
        current_pos[1] += self.r*np.sin(theta)*np.sin(phi)
        current_pos[2] += self.r*np.cos(theta)
        # Update list of visited points
        self.visited_points.append(tuple(current_pos))

    def walk_with_self_avoid(self,nsteps=100):
        """Walk nsteps steps of self-avoiding random walk"""
        for i in range(nsteps):
            self.walk_one_step()
            print(self.visited_points)
            current_pos = np.array(self.visited_points[-1])
            visited_points = np.array(self.visited_points[:-1]) # This may be very ineffective programming
            if any(np.sqrt(sum((t-current_pos)**2)) < 2*self.R for t in visited_points):
                print('Managed to walk',len(visited_points)-1,'steps')
                break

# gridwalk = Grid_walker()
# # print(gridwalk.get_multiple_end_to_end_distances(nwalks=10))
# gridwalk.walk_without_avoid()
# # gridwalk.walk_with_self_avoid(nsteps=10)
# print(gridwalk.visited_points)
# gridwalk.plot_the_walk()

chainwalk = Freely_jointed_chain()
chainwalk.walk_without_avoid(nsteps=10)
# chainwalk.walk_with_self_avoid()
print(chainwalk.visited_points)
chainwalk.plot_the_walk()
