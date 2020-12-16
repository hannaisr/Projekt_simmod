import random as rnd
import numpy as np
import matplotlib.pyplot as plt # For 3D plot
from mpl_toolkits.mplot3d import Axes3D # For 3D plot

class Walker():
    """Walked positions are stored in a list"""
    # Initiate list of visited points
    origin = [0,0,0]
    visited_points = [origin]

    # Define step length
    r = 1

    # Define bead radius, at largest 1/2*r
    rho = 0.2

    # Store last walking direction
    last_direction = 0
    last_directionList = [0,0,0]

    def walk_one_step(self, limited=False):
        """Implemented by child class"""
        pass

    def walk_without_avoid(self,nsteps=100, limited=False):
        """Walk nsteps steps of non-self-avoiding random walk."""
        self.restart()
        for i in range(nsteps):
            self.walk_one_step(limited)

    def walk_with_self_avoid(self,nsteps=100,limited=True,cheat=False):
        """Walk nsteps steps of self-avoiding random walk. If cheat=True, each step is repeated until it is successful or it has failed 100 times."""
        if cheat==True:
            tries_per_step = 100 # Maximum number of times to try again if the step is unacceptable
        else:
            tries_per_step = 1

        self.restart()
        try_again = True
        # Try to assemble a sequence until successful
        while try_again is True:
            try_again = False
            for i in range(nsteps):
                for j in range(tries_per_step):
                    self.walk_one_step(limited)
                    # Test if the site is already occupied.
                    try_again=self.test_avoid()
                    if try_again is True:   # Remove last site and start again
                        self.visited_points.pop()
                    else:
                        break
                # In case of self interception, break attempt immediately
                if try_again is True:
                    print('Managed to walk',len(self.visited_points)-2,'steps')
                    self.restart()
                    break
        print('Managed to walk', len(self.visited_points) - 1, 'steps')

    def test_avoid(self):
        """Implemented by child class"""
        pass

    def restart(self):
        """Resets list of visited points."""
        self.visited_points = [self.origin]
        self.last_direction = 0

    def plot_the_walk(self,beads=False):
        """Plots the walk in 3D."""
        x = [i[0] for i in self.visited_points]
        y = [i[1] for i in self.visited_points]
        z = [i[2] for i in self.visited_points]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x,y,z)
        if beads is True:
            cmap = get_cmap(len(x))
            for i in range(len(x)):
                phi, theta = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
                x_sphere = x[i] + self.rho * np.cos(phi) * np.sin(theta)
                y_sphere = y[i] + self.rho * np.sin(phi) * np.sin(theta)
                z_sphere = z[i] + self.rho * np.cos(theta)
                ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color=cmap(i))
        plt.show()

    def get_end_to_end_distance(self):
        """Calculate end-to-end distance of already walked walk."""
        last_point = self.visited_points[-1]
        return np.sqrt((last_point[0]-self.origin[0])**2+(last_point[1]-self.origin[1])**2+(last_point[2]-self.origin[2])**2)

    def get_multiple_end_to_end_distances(self,nsteps=100,nwalks=10,avoid=False):
        """Returns a list of end-to-end distances for nwalks number of walks of length nsteps"""
        etedist_list = np.zeros(nwalks)
        for i in range(nwalks):
            self.restart()
            if avoid is True:
                self.walk_with_self_avoid(nsteps=nsteps)
            else:
                self.walk_without_avoid(nsteps=nsteps)
            etedist_list[i] = self.get_end_to_end_distance()
        return etedist_list

    def plot_multiple_end_to_end_distances(self,nwalks=10,avoid=False):
        """Plots end-to-end distance RMS, RMS fluctuation and standard error estimate for nwalks walks by number of steps"""
        rms=[]
        rms_fluc = []
        std_err = []
        step_numbers = range(100,1100,100)
        for nsteps in step_numbers:
            etedist_list = self.get_multiple_end_to_end_distances(nsteps=nsteps,nwalks=nwalks,avoid=avoid)
            #RMS end-to-end distance
            rms.append(np.sqrt(np.mean(np.square(etedist_list))))
            #RMS fluctuation estimate
            rms_fluc.append(np.sqrt((np.mean(np.square(etedist_list))-np.mean(etedist_list)**2)*nwalks/(nwalks-1)))
            #Standard error estimate
            std_err.append(np.sqrt((np.mean(np.square(etedist_list))-np.mean(etedist_list)**2)/(nwalks-1)))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(step_numbers,rms,label="RMS End-to-End Distance")
        ax.plot(step_numbers,rms_fluc,label="RMS Fluctuation Estimate")
        ax.plot(step_numbers,std_err,label="Standard Error Estimate")
        ax.set_xlabel("Number of steps")
        plt.legend()
        plt.show()

class Grid_walker(Walker):
    def walk_one_step(self, limited=False):
        possible_directions = [-3,-2,-1,1,2,3]
        current_pos = self.visited_points[-1][:]
        # Get walking direction
        direction = rnd.choice(possible_directions)
        if limited == True:
            while direction == -self.last_direction:
                direction = rnd.choice(possible_directions)
        self.last_direction = direction
        # Update the coordinates
        if direction == 1:
            current_pos[0] += self.r
        elif direction == -1:
            current_pos[0] -= self.r
        elif direction == 2:
            current_pos[1] += self.r
        elif direction == -2:
            current_pos[1] -= self.r
        elif direction == 3:
            current_pos[2] += self.r
        elif direction == -3:
            current_pos[2] -= self.r
        # Update list of visited points
        self.visited_points.append(current_pos)

    def test_avoid(self):
        """Test if latest site is already occupied. Return True if so, False if not."""
        if any(t == self.visited_points[-1] for t in self.visited_points[:-1]):
            return True
        return False

class Freely_jointed_chain(Walker):
    def walk_one_step(self, limited=False):
        current_pos = self.visited_points[-1][:]
        # Get walking direction
        if limited == True and self.last_direction != 0:
            # Define direction to walk back the same way
            theta_back = np.pi-self.last_direction[0]
            phi_back = np.pi+self.last_direction[1]
            # Define angle for resctriction cone
            alpha = 2*np.arsin(self.rho/self.r)
            # Define accepted angles in order to not walk back the same way
            accepted_theta = [theta_back+alpha,2*np.pi+theta_back-alpha]
            accepted_phi =  [phi_back+alpha,2*np.pi+phi_back-alpha]
            # Get new direction
            theta = rnd.uniform(accepted_theta[0],accepted_theta[1])
            phi = rnd.uniform(accepted_phi[0],accepted_phi[1])
        else:
            theta = rnd.uniform(0,np.pi)
            phi = rnd.uniform(0,2*np.pi)
        self.last_direction = [theta,phi]
        # Update the coordinates
        current_pos[0] += self.r*np.sin(theta)*np.cos(phi)
        current_pos[1] += self.r*np.sin(theta)*np.sin(phi)
        current_pos[2] += self.r*np.cos(theta)
        # Update list of visited points
        self.visited_points.append(current_pos)

    def test_avoid(self):
        """Test if latest site is already occupied - continuous case"""
        #The distance between successive sphere centres is self.r. Interception between any two spheres occurs if their centres are less apart than their diameter
        for point in self.visited_points[:-1]:
            r_centres = np.sqrt((point[0] - self.visited_points[-1][0])**2 + (point[1] - self.visited_points[-1][1])**2 + (point[2] - self.visited_points[-1][2])**2)
            if r_centres < 2*self.rho:
                # Self-intercept - needs to restart the process
                return True
        return False

class Directed_walker(Freely_jointed_chain):
    """Walks in specific direction with specified distribution."""
    def walk_one_step(self, limited=False):
        current_pos = self.visited_points[-1][:]
        # Get walking direction
        theta = rnd.normalvariate(np.pi/2,np.pi/4)
        phi = rnd.normalvariate(np.pi,np.pi/2)
        self.last_direction = [theta,phi]
        # Update the coordinates
        current_pos[0] += self.r*np.sin(theta)*np.cos(phi)
        current_pos[1] += self.r*np.sin(theta)*np.sin(phi)
        current_pos[2] += self.r*np.cos(theta)
        # Update list of visited points
        self.visited_points.append(current_pos)


class Grid_walker_stepl_variations(Grid_walker, Freely_jointed_chain):
    def walk_one_step(self, limited=False):
        self.r = rnd.normalvariate(1,1) # Varation in step length
        Grid_walker.walk_one_step(self,limited)

    def test_avoid(self):
        Freely_jointed_chain.test_avoid(self)

class Brownian(Walker):
    """3D chain generated by Brownian motion in each spatial direction"""
    def walk_one_step(self, limited=False):
        current_pos = self.visited_points[-1][:]
        # Get walking direction
        direction = [rnd.normalvariate(0,1) for i in range(3)] #TODO: How to make it have expected step length of 1?
        if limited == True:
            while direction == -1*self.last_directionList: #TODO: Improve to follow the Freely Jointed Chain method
                direction = [rnd.normalvariate(1/np.sqrt(3),1) for i in range(3)]
        self.last_direction = direction
        # Update the coordinates
        current_pos[0] += direction[0]
        current_pos[1] += direction[1]
        current_pos[2] += direction[2]
        # Update list of visited points
        self.visited_points.append(current_pos)

    def test_avoid(self):
        """Test if latest site is already occupied - continuous case"""
        #Interception between any two spheres occurs if their centres are less apart than their diameter
        for point in self.visited_points[:-1]:
            r_centres = np.sqrt((point[0] - self.visited_points[-1][0])**2 + (point[1] - self.visited_points[-1][1])**2 + (point[2] - self.visited_points[-1][2])**2)
            if r_centres < 2*self.rho:
                # Self-intercept - needs to restart the process
                return True
        return False

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n + 1)

def main():
    # gridwalk = Grid_walker()
    # gridwalk.plot_multiple_end_to_end_distances()
    # print(gridwalk.get_multiple_end_to_end_distances(nwalks=10,avoid=False))
    # gridwalk.walk_without_avoid(nsteps=100,limited=False)
    # gridwalk.walk_with_self_avoid(nsteps=100,limited=True)
    # gridwalk.plot_the_walk(beads=False)

    # chainwalk = Freely_jointed_chain()
    # chainwalk.plot_multiple_end_to_end_distances(nwalks=100)
    # chainwalk.walk_without_avoid(nsteps=1000,limited=False)
    # chainwalk.plot_the_walk(beads=False)
    # chainwalk.walk_with_self_avoid(nsteps=10000,limited=True)
    # chainwalk.plot_the_walk(beads=False)

    # dirwalk = Directed_walker()
    # dirwalk.walk_without_avoid(nsteps=1000)
    # dirwalk.plot_the_walk(beads=False)
    # dirwalk.plot_multiple_end_to_end_distances(avoid=False)

    # brownian = Brownian()
    # brownian.walk_without_avoid(nsteps=300)
    # brownian.walk_with_self_avoid(nsteps=50)
    # brownian.plot_the_walk(beads=True)
    # brownian.plot_multiple_end_to_end_distances(avoid=False)

    grid_walker_stepl_variations = Grid_walker_stepl_variations()
    grid_walker_stepl_variations.walk_without_avoid(nsteps=10)
    grid_walker_stepl_variations.plot_the_walk(beads=False)

main()
