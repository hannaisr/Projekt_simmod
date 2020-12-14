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

    def walk_one_step(self):
        """Implemented by child class"""
        pass

    def walk_without_avoid(self,nsteps=100, limited=False):
        """Walk nsteps steps of non-self-avoiding random walk."""
        self.restart()
        for i in range(nsteps):
            self.walk_one_step(limited)

    def walk_with_self_avoid(self,nsteps=100):
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
            current_pos[0] += 1
        elif direction == -1:
            current_pos[0] -= 1
        elif direction == 2:
            current_pos[1] += 1
        elif direction == -2:
            current_pos[1] -= 1
        elif direction == 3:
            current_pos[2] += 1
        elif direction == -3:
            current_pos[2] -= 1
        # Update list of visited points
        self.visited_points.append(current_pos)

    def walk_with_self_avoid(self,nsteps=100,limited=True):
        """Walk nsteps steps of self-avoiding random walk"""
        self.restart()
        try_again = True
        while try_again is True:
            try_again = False
            for i in range(nsteps):
                self.walk_one_step(limited)
                #In case of self-interception, abort attempt and retry
                if self.visited_points[-1] in self.visited_points[:-1]:
                    print('Managed to walk',len(self.visited_points)-1,'steps')
                    self.restart()
                    try_again = True
                    break
        print('Managed to walk', len(self.visited_points) - 1, 'steps')

class Freely_jointed_chain(Walker):
    def walk_one_step(self, limited=False):
        current_pos = self.visited_points[-1][:]
        # Get walking direction
        if limited == True and self.last_direction != 0:
            # Define direction to walk back the same way
            theta_back = np.pi-self.last_direction[0]
            phi_back = np.pi+self.last_direction[1]
            # Define accepted angles in order to not walk back the same way
            accepted_theta = [theta_back+self.rho,2*np.pi+theta_back-self.rho]  # - TODO - intervallet är nu för stort och kan göras mer exakt.
            accepted_phi =  [phi_back+self.rho,2*np.pi+phi_back-self.rho]
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

    def walk_with_self_avoid(self,nsteps=100,limited=True):
        """Walk nsteps steps of self-avoiding random walk"""
        self.restart()
        try_again = True
        # Try to assemble a sequence until successful
        while try_again is True:
            try_again = False
            for i in range(nsteps):
                self.walk_one_step(limited)
                # Test if the site is already occupied
                try_again = self.test_avoid()
                # In case of self interception, break attempt immediately
                if try_again is True:
                    print('Managed to walk',len(self.visited_points)-1,'steps')
                    self.restart()
                    break
        print('Managed to walk', len(self.visited_points) - 1, 'steps')

    def test_avoid(self):
        """Test if latest site is already occupied"""
        #The distance between neighboring sphere centres is self.r, so each sphere has radius 1/2*self.r
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
    gridwalk = Grid_walker()
    # gridwalk.plot_multiple_end_to_end_distances()
    # print(gridwalk.get_multiple_end_to_end_distances(nwalks=10,avoid=False))
    # gridwalk.walk_without_avoid(nsteps=100,limited=True)
    # gridwalk.walk_with_self_avoid(nsteps=50)
    # gridwalk.plot_the_walk(beads=False)

    chainwalk = Freely_jointed_chain()
    # chainwalk.plot_multiple_end_to_end_distances(nwalks=100)
    # chainwalk.walk_without_avoid(nsteps=100,limited=True)
    chainwalk.walk_with_self_avoid(nsteps=100,limited=True)
    chainwalk.plot_the_walk(beads=False)

main()
