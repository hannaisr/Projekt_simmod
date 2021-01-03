import random as rnd
import numpy as np
import matplotlib.pyplot as plt # For 3D plot
from mpl_toolkits.mplot3d import Axes3D # For 3D plot
from scipy.stats import norm # histogram fitting
from scipy import stats
import time

class Walker():
    """Walked positions are stored in a list"""
    # Initiate list of visited points
    rho0 = 0.499 # size of first bead, the very least 1/2 step length.
    # Modify generate_rho() to manage method for generating the sizes of the other beads
    origin = [0,0,0,rho0] # position and bead size are stored in list
    visited_points = [origin]
    length = 0
    my=1
    sigma=0.1

    # Define step length
    r = my

    # Store last walking direction
    last_direction = 0

    variate_rho=False

    def walk_one_step(self, limited=False):
        """Implemented by child class"""
        pass

    def walk_without_avoid(self,nsteps=100, limited=False):
        """Walk nsteps steps of non-self-avoiding random walk."""
        self.restart()
        for i in range(nsteps):
            self.walk_one_step(limited)

    def walk_with_self_avoid_forced(self,nsteps=100,limited=True):
        """Walk nsteps steps of self-avoiding random walk, redoing each step until it is successful or it has failed tries_per_step times."""

        tries_per_step = 100 # Maximum number of times to try again if the step is unacceptable

        self.restart()
        try_again = True
        # Try to assemble a sequence until successful
        while try_again is True:
            try_again = False
            for i in range(nsteps):
                for j in range(tries_per_step):
                    self.walk_one_step(limited)
                    # Test if the site is already occupied, or if the step length is less than the bead diameter.
                    try_again=self.test_avoid()
                    if try_again is True:   # Remove last site and start again
                        self.visited_points.pop()
                    else:
                        break
                # In case of self interception, break attempt immediately
                if try_again is True:
                    # print('Managed to walk',len(self.visited_points)-2,'steps')
                    self.restart()
                    break
        # print('Managed to walk', len(self.visited_points) - 1, 'steps')

    def walk_with_self_avoid(self,nsteps=100,limited=True):
        """Walk nsteps steps of self-avoiding random walk. Returns the number of walks that failed."""

        nfails = 0

        self.restart()
        try_again = True
        # Try to assemble a sequence until successful
        while try_again is True:
            try_again = False
            for i in range(nsteps):
                self.walk_one_step(limited)
                # Test if the site is already occupied.
                try_again=self.test_avoid()
                # In case of self interception, break attempt immediately
                if try_again is True:
                    # print('Managed to walk',len(self.visited_points)-2,'steps')
                    nfails += 1
                    self.restart()
                    break
        # print('Managed to walk', len(self.visited_points) - 1, 'steps')
        return nfails

    def generate_rho(self):
        if self.variate_rho is True:
            return(self.r-self.visited_points[-1][3]-0.0000001) # Current and last bead radiuses alwys cover the whole length of the step
        else:
            return(self.rho0)

    def success_rate(self,nsteps=20,limited=True):
        """Calculates success rate from nsuccessful_walks number of successful walks."""
        total_fails = 0
        nsuccessful_walks = 100
        for i in range(nsuccessful_walks):
            total_fails += self.walk_with_self_avoid(nsteps,limited)
        return nsuccessful_walks/(nsuccessful_walks+total_fails)


    def test_avoid(self):
        """Implemented by child class"""
        pass

    def restart(self):
        """Resets list of visited points."""
        self.visited_points = [self.origin]
        self.last_direction = 0
        self.length = 0
        self.r = self.my

    def plot_the_walk(self,beads=False):
        """Plots the walk in 3D."""
        x = [i[0] for i in self.visited_points]
        y = [i[1] for i in self.visited_points]
        z = [i[2] for i in self.visited_points]
        rho = [i[3] for i in self.visited_points]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x,y,z)
        if beads is True:
            cmap = get_cmap(len(x))
            for i in range(len(x)):
                phi, theta = np.mgrid[0:2 * np.pi:200j, 0:np.pi:100j]
                x_sphere = x[i] + rho[i] * np.cos(phi) * np.sin(theta)
                y_sphere = y[i] + rho[i] * np.sin(phi) * np.sin(theta)
                z_sphere = z[i] + rho[i] * np.cos(theta)
                ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color=cmap(i))
        plt.show()

    def get_end_to_end_distance(self):
        """Calculate end-to-end distance of already walked walk."""
        last_point = self.visited_points[-1]
        return np.sqrt((last_point[0]-self.origin[0])**2+(last_point[1]-self.origin[1])**2+(last_point[2]-self.origin[2])**2)

    def get_multiple_end_to_end_distances(self,nsteps=100,nwalks=10,avoid=False,limited=True,forced=False):
        """Returns a list of end-to-end distances and chain lengths for nwalks number of walks of length nsteps"""
        etedist_list = np.zeros(nwalks)
        length_list = np.zeros(nwalks)
        for i in range(nwalks):
            self.restart()
            if avoid is True:
                if forced is True:
                    self.walk_with_self_avoid_forced(nsteps=nsteps,limited=limited)
                else:
                    self.walk_with_self_avoid(nsteps=nsteps,limited=limited)
            else:
                self.walk_without_avoid(nsteps=nsteps)
            etedist_list[i] = self.get_end_to_end_distance()
            length_list[i] = self.length
            # print("Finished",i+1,"walks")
        return etedist_list, length_list

    def plot_multiple_end_to_end_distances(self,nwalks=10,avoid=False,limited=False,forced=False,holdon=False):
        """Plots end-to-end distance RMS, RMS fluctuation and standard error estimate for nwalks walks by chain length"""
        rms=[]
        rms_fluc = []
        std_err = []
        step_numbers = range(5,14,1)
        for nsteps in step_numbers:
            print(nsteps)
            etedist_list, length_list = self.get_multiple_end_to_end_distances(nsteps=nsteps,nwalks=nwalks,avoid=avoid,limited=limited,forced=forced)
            #RMS end-to-end distance
            rms.append(np.sqrt(np.mean(np.square(etedist_list))))
            #RMS fluctuation estimate
            rms_fluc.append(np.sqrt((np.mean(np.square(etedist_list))-np.mean(etedist_list)**2)*nwalks/(nwalks-1)))
            #Standard error estimate
            std_err.append(np.sqrt((np.mean(np.square(etedist_list))-np.mean(etedist_list)**2)/(nwalks-1)))
        chain_lengths = [i * self.r for i in step_numbers]
        if holdon==True:
            return chain_lengths,rms,rms_fluc,std_err
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(chain_lengths,rms,label="RMS End-to-End Distance")
            ax.plot(chain_lengths,rms_fluc,label="RMS Fluctuation Estimate")
            ax.plot(chain_lengths,std_err,label="Standard Error Estimate")
            ax.set_xlabel("Chain length")
            plt.legend()

            ext=""
            if avoid==True:
                ext+=", self-avoiding"
            if limited==True:
                ext += ", limited"
            if forced == True:
                ext+=", forced"

            plt.suptitle("End-to-end distance measures vs chain length \n for " +str(nwalks)+" walks of "+self.name+ext)
            plt.show()

    def plot_multiple_end_to_end_distances_holdon(self,nwalks=10):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel("Chain length")
        plt.suptitle("End-to-end distance measures vs chain length \n for " +str(nwalks)+" walks of "+self.name+', bead radius='+str(self.rho0))
        cmap=get_cmap(n=4)

        #Non-self avoid
        i=0
        chain_lengths, rms, rms_fluc, std_err = self.plot_multiple_end_to_end_distances(nwalks=nwalks,avoid=False,limited=False,forced=False,holdon=True)
        ax.plot(chain_lengths, rms, label="Non self avoiding",color=cmap(i))
        ax.plot(chain_lengths, rms_fluc,color=cmap(i))
        ax.plot(chain_lengths, std_err,color=cmap(i))

        #Self avoid
        i=1
        chain_lengths, rms, rms_fluc, std_err = self.plot_multiple_end_to_end_distances(nwalks=nwalks,avoid=True,limited=False,forced=False,holdon=True)
        ax.plot(chain_lengths, rms, label="Self avoiding",color=cmap(i))
        ax.plot(chain_lengths, rms_fluc,color=cmap(i))
        ax.plot(chain_lengths, std_err,color=cmap(i))

        #Self avoid:limited
        i=2
        chain_lengths, rms, rms_fluc, std_err = self.plot_multiple_end_to_end_distances(nwalks=nwalks,avoid=True,limited=True,forced=False,holdon=True)
        ax.plot(chain_lengths, rms, label="Self avoiding: limited",color=cmap(i))
        ax.plot(chain_lengths, rms_fluc,color=cmap(i))
        ax.plot(chain_lengths, std_err,color=cmap(i))

        #Self avoid:forced
        i=3
        chain_lengths, rms, rms_fluc, std_err = self.plot_multiple_end_to_end_distances(nwalks=nwalks,avoid=True,limited=True,forced=True,holdon=True)
        ax.plot(chain_lengths, rms, label="Self avoiding: forced",color=cmap(i))
        ax.plot(chain_lengths, rms_fluc,color=cmap(i))
        ax.plot(chain_lengths, std_err,color=cmap(i))


        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        textstr="RMS"
        # place a text box in upper left in axes coords
        ax.text(0.90, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
        textstr = "RMS\nfluctuation"
        # place a text box in upper left in axes coords
        ax.text(0.90, 0.3, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
        textstr = "Standard error\nestimate"
        # place a text box in upper left in axes coords
        ax.text(0.90, 0.1, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)

        plt.legend()
        plt.show()

    def hist_quotient_length_etedist(self,nsteps=100,nwalks=10,avoid=False,limited=True,forced=False):
        """Plots the quotient between total chain length and end-to-end distance"""
        etedist_list, length_list = self.get_multiple_end_to_end_distances(nsteps=nsteps,nwalks=nwalks,avoid=avoid,limited=limited,forced=forced)
        quotients = [etedist_list[i]/length_list[i] for i in range(len(length_list))]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(x=quotients,bins=15)
        ax.vlines(ymin=0,ymax=nwalks/5,x=np.mean(quotients),color='red',label="Mean")
        ax.set_xlabel("End-to-End distance/length")
        plt.legend()
        plt.show()

    def plot_etedist_normal(self,nwalks=1000,avoid=False,limited=False,forced=False):
        """Plots fitted mu and sigma for end-to-end distance vs length of chain."""
        etedist_lists = []
        nsteps_list = np.arange(0,30,1)
        length_list = list(nsteps_list*self.my)
        mu_list = []
        sigma_list = []
        # Get end-to-end distances
        for nsteps in nsteps_list:
            etedist_list, _ = self.get_multiple_end_to_end_distances(nsteps=nsteps,nwalks=nwalks,avoid=avoid,limited=limited,forced=forced)
            [mu,sigma] = norm.fit(etedist_list)
            mu_list.append(mu)
            sigma_list.append(sigma)
            print(nsteps,"nsteps finished")
        if limited:
            title = "Normal distribution parameters vs mean chain length\n "+str(self.name)+", "+str(nwalks)+" walks per chain length, limited,\n "
        else:
            title = "Normal distribution parameters vs mean chain length\n "+str(self.name)+", "+str(nwalks)+" walks per chain length, not limited,\n "
        if forced is True:
            title += "forced "
        if avoid:
            title += "self avoiding"
        else:
            title += "not self avoiding"
        labels = ["mu","sigma"]
        plot2D(title,"Chain length",None,[length_list,length_list],[mu_list,sigma_list],labels)
        print("List of mus:",mu_list,"\n List of sigmas:",sigma_list)
        return mu_list, sigma_list, length_list

    def plot_etedist_parameters_multiple_versions(self):
        """Plots end-to-end distance mu and sigma for different methods of random walk"""
        mu1,sigma1,lenlist1 = self.plot_etedist_normal(limited=False)
        mu3,sigma3,lenlist3 = self.plot_etedist_normal(avoid=True,limited=True)
        mu4,sigma4,lenlist4 = self.plot_etedist_normal(limited=False,forced=True,avoid=True)
        mu2,sigma2,lenlist2 = self.plot_etedist_normal(limited=False,avoid=True,forced=False)
        x = range(0,30,1)
        plot2D("End-to-end distance mu vs chain lenght\n for "+str(self.name),"Chain length", "End-to-end distance mu", [lenlist1,lenlist2,lenlist3,lenlist4],[mu1,mu2,mu3,mu4],["Not limited, not avoiding","Self avoiding", "Self avoiding limited","Forced self avoiding"])
        plot2D("End-to-end distance sigma vs chain lenght\n for "+str(self.name),"Chain length", "End-to-end distance sigma", [lenlist1,lenlist2,lenlist3,lenlist4],[sigma1,sigma2,sigma3,sigma4],["Not limited, not avoiding","Self avoiding", "Self avoiding limited","Forced self avoiding"])
        plt.show()

    def plot_success_rate_vs_nsteps(self,step_numbers=range(1,25,5),limited=True):
        """Plots success rate to number of steps in walk"""
        #step_numbers = range(10,120,10)   # Grid
        step_numbers = range(1,25,5)        # Freely jointed
        success_rates = []
        for nsteps in step_numbers:
            success_rates.append(self.success_rate(nsteps,limited))
            print(success_rates)
        title = "Success rate vs number of steps for \n"+str(self.name)+", "
        if limited is True:
            title += 'limited'
        else:
            title += 'not limited'
        plot2D(title,'Number of steps','Success rate',step_numbers,success_rates)
        plt.show()
        return step_numbers,success_rates

    def plot_bead_size_variation(self,nsteps=100,nwalks=10,my=True,success=False,limited=False,forced=False):
        """Plots the relationship between quotient bead size/my or sigma of step length and RMS End-to-End distance or success rate in self-avoiding chains"""
        qs=[] #Quotients bead size/expected value of step length
        rms=[]
        rms_fluc = []
        std_err = []
        success_rates=[]
        ext = ""
        if limited == True:
            ext += ", limited"
        if forced == True:
            ext += ", forced"
        rhos = list(range(1,5))
        rhos.append(4.99)
        for rho in rhos:
            print(rho)
            self.rho0=rho/10
            if my==True: #Bead size by expected value
                qs.append(2*self.rho0/self.my)
            else: #Bead size by standard deviation
                qs.append(2*self.rho0/self.sigma)
            if success==False:#y=RMS
                etedist_list, length_list = self.get_multiple_end_to_end_distances(nsteps=nsteps, nwalks=nwalks, avoid=True, limited=limited,forced=forced)
                # RMS end-to-end distance
                rms.append(np.sqrt(np.mean(np.square(etedist_list))))
                # RMS fluctuation estimate
                rms_fluc.append(np.sqrt((np.mean(np.square(etedist_list)) - np.mean(etedist_list) ** 2) * nwalks / (nwalks - 1)))
                # Standard error estimate
                std_err.append(np.sqrt((np.mean(np.square(etedist_list)) - np.mean(etedist_list) ** 2) / (nwalks - 1)))
            else:#y=Success rate
                success_rates.append(self.success_rate(nsteps=nsteps,limited=limited))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if my==True:
            ax.set_xlabel("Bead diameter/Expected value of step length")
        else:
            ax.set_xlabel("Bead diameter/Standard deviation of step length")
        if success == False:
            ax.plot(qs, rms, label="RMS End-to-End Distance")
            ax.plot(qs, rms_fluc, label="RMS Fluctuation Estimate")
            ax.plot(qs, std_err, label="Standard Error Estimate")
            plt.suptitle("End-to-end distance measures vs bead size \n for " + str(nwalks) +" walks of "+ str(nsteps) + "-steps "+ self.name+ext)
        else:
            ax.plot(qs,success_rates)
            ax.set_ylabel("Success rate")
            plt.suptitle("Success rate vs bead size for " + str(nsteps) + "-steps \n"+self.name+ext)
        plt.legend()
        plt.show()

    def normplot(self,nsteps=100,nwalks=10,avoid=False,limited=True,forced=False):
        etedistlist,l = self.get_multiple_end_to_end_distances(nsteps=nsteps,nwalks=nwalks,avoid=avoid,limited=limited,forced=forced)

        # Calculate quantiles and least-square-fit curve
        (quantiles, values), (slope, intercept, r) = stats.probplot(etedistlist, dist='norm')

        # plot results
        plt.plot(values, quantiles, 'ob')
        plt.plot(quantiles * slope + intercept, quantiles, 'r')

        # define ticks
        ticks_perc = [1, 5, 10, 20, 50, 80, 90, 95, 99]

        # transfrom them from precentile to cumulative density
        ticks_quan = [stats.norm.ppf(i / 100.) for i in ticks_perc]

        # assign new ticks
        plt.yticks(ticks_quan, ticks_perc)

        # show plot
        plt.grid()
        plt.show()

class Grid_walker(Walker):
    def __init__(self):
        self.name = "Grid walker, r="+str(self.r)

    def walk_one_step(self, limited=False):
        possible_directions = [-3,-2,-1,1,2,3]
        current_pos = self.visited_points[-1][:]
        current_pos[3] = self.generate_rho()
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
        self.length += self.r

    def test_avoid(self):
        """Test if latest site is already occupied. Return True if so, False if not."""
        # if self.r < 2*self.rho:
        #     return True
        if self.visited_points[0][:3] == self.visited_points[-1][:3]:
            return True
        elif any(t == self.visited_points[-1][:3] for t in [i[:3] for i in self. visited_points[:-1]]):
            return True
        return False

class Freely_jointed_chain(Walker):
    def __init__(self):
        self.name = 'Freely jointed chain, r='+str(self.r)

    def walk_one_step(self, limited=False):
        current_pos = self.visited_points[-1][:]
        # Get walking direction and bead size
        theta = rnd.uniform(0,np.pi)
        phi = rnd.uniform(0,2*np.pi)
        rho = self.generate_rho()
        if limited == True and self.last_direction != 0:
            # Define direction to walk back the same way
            theta_back = np.pi-self.last_direction[0]
            phi_back = np.pi+self.last_direction[1]
            while (self.r*np.sin(theta)*np.cos(phi)-self.r*np.sin(theta_back)*np.cos(phi_back))**2+(self.r*np.sin(theta)*np.sin(phi)-self.r*np.sin(theta_back)*np.sin(phi_back))**2+(self.r*np.cos(theta)-self.r*np.cos(theta_back))**2 < (current_pos[3]+rho)**2:
                theta = rnd.uniform(0,np.pi)
                phi = rnd.uniform(0,2*np.pi)
        self.last_direction = [theta,phi]
        # Update the coordinates
        current_pos[0] += self.r*np.sin(theta)*np.cos(phi)  # x
        current_pos[1] += self.r*np.sin(theta)*np.sin(phi)  # y
        current_pos[2] += self.r*np.cos(theta)              # z
        current_pos[3] = rho
        # Update list of visited points
        self.visited_points.append(current_pos)
        self.length += self.r

    def test_avoid(self):
        """Test if latest site is already occupied - continuous case"""
        cp = self.visited_points[-1]    # Current position
        #The distance between successive sphere centres is self.r. Interception between any two spheres occurs if their centres are less apart than their diameter
        # if self.r < 2*self.rho:
        #     return True
        for point in self.visited_points[:-1]:
            r_centres = np.sqrt((point[0] - cp[0])**2 + (point[1] - cp[1])**2 + (point[2] - cp[2])**2)
            if r_centres < point[3]+cp[3]:
                # Self-intercept - needs to restart the process
                return True
        return False
        # TODO Implement method for avoiding previous _paths_, not just previous positions.

class Reptation_walker(Grid_walker):
    origin = [0,0,0]

    """Reptation algorithm to generate SAWs. Algoritm on hiskp.uni-bonn... pg. 4"""
    def __init__(self,nsteps=100,name='Reptation walker'):
        # Generate a basis self-avoiding walk
        self.gridwalk = Grid_walker()
        self.gridwalk.walk_with_self_avoid(nsteps=nsteps)
        self.visited_points = self.gridwalk.visited_points #The gridwalk.visited_points will be modified as saw is
        self.nsteps=nsteps
        self.name = name

    def test_avoid(self):
        """Test if latest site is already occupied. Return True if so, False if not."""
        if self.visited_points[0][:3] == self.visited_points[-1][:3]:
            return True
        elif any(t == self.visited_points[-1][:3] or t == self.visited_points[0][:3] for t in [i[:3] for i in self.visited_points[1:-1]]):
            return True
        return False

    def walk_with_self_avoid(self,nsteps=100,limited=True):
        """Generates a new SAW using the reptation model based on the initialized SAW"""
        #TODO: The signatures of the methods aren't optimal since the methods are slightly different, but it works
        #Save a copy of the old saw #TODO: Performance improvement possible if computational cost is an issue
        old_visited_points = list(self.visited_points)
        #Length nsteps defined in class attribute
        self.length = self.nsteps*self.r

        #***Algorithm start***
        #Choose an end point at random
        choice = rnd.choice([0, -1])
        #Remove this end
        self.visited_points.pop(choice) #One end
        #Add a step on the other end identified as "current_pos"
        current_pos = list(self.visited_points[choice * (-1) - 1])  # The other end
        prev_pos = self.visited_points[choice * (-1) + (choice + 1) * (-2)] #om choice i -1: prev_pos index= 1; om choice i 0: prev_pos index= -2
        # Get walking direction
        possible_directions = [-3, -2, -1, 1, 2, 3]
        direction = rnd.choice(possible_directions)
        #Get the last direction
        for i in range(3):
            step = current_pos[i] - prev_pos[i] #Step is either +1 or -1
            if abs(step) > 0:
                last_direction = (i + 1) * step
                break
        while direction == -last_direction:
            direction = rnd.choice(possible_directions)
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
        if choice == -1:
            self.visited_points.insert(0, current_pos)
        else:
            self.visited_points.append(current_pos)
        # In case of self-intersection: revert to previous configuration. Else, retain new configuration #TODO: Is this the right interpretation?
        if self.test_avoid() == True:
            self.visited_points = old_visited_points
        self.gridwalk.visited_points = self.visited_points
        return 0

    def get_multiple_end_to_end_distances(self,nsteps=100,nwalks=10,avoid=False,limited=True):
        """Returns a list of end-to-end distances and chain lengths for nwalks number of walks"""
        etedist_list = np.zeros(nwalks)
        length_list = np.zeros(nwalks)
        for i in range(nwalks):
            self.walk_with_self_avoid()
            etedist_list[i] = self.get_end_to_end_distance()
            length_list[i] = self.length
        return etedist_list, length_list

class Directed_walker(Freely_jointed_chain):
    """Walks in specific direction with specified distribution."""
    def __init__(self,name='Freely jointed chain with directed walk'):
        self.name=name

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
        self.length += self.r

class Grid_walker_stepl_variations(Grid_walker, Freely_jointed_chain):
    """Parallell/perpendicular motion with a randomly distributed step length"""
    def __init__(self,distribution="N"):
        self.distribution = distribution
        self.name = "Grid walker with step length variation, mu="+str(self.my)+", sigma="+str(self.sigma)

    def walk_one_step(self, limited=False):
        if self.distribution == "N":
            self.r = rnd.normalvariate(self.my,self.sigma) # Varation in step length
        elif self.distribution == "exp":
            self.r = rnd.expovariate(1/self.my) # Varation in step length
        Grid_walker.walk_one_step(self,limited)

    def test_avoid(self):
        return Freely_jointed_chain.test_avoid(self)

class Freely_jointed_chain_stepl_variations(Freely_jointed_chain):
    """Freely jointed chain with randomly distributed step lengths"""
    def __init__(self, distribution="N"):
        self.distribution = distribution
        self.name = "Freely jointed chain with step length variations, mu="+str(self.my)+", sigma="+str(self.sigma)

    def walk_one_step(self, limited=False):
        if self.distribution == "N":
            self.r = rnd.normalvariate(self.my,self.sigma)  # Varation in step length
        elif self.distribution == "exp":
            self.r = rnd.expovariate(1/self.my)  # Varation in step length
        Freely_jointed_chain.walk_one_step(self,limited)

def reptation_plot_multiple_end_to_end_distances(nwalks=10,avoid=False,limited=False):
    """Plots end-to-end distance RMS, RMS fluctuation and standard error estimate for nwalks walks by chain length"""
    rms=[]
    rms_fluc = []
    std_err = []
    step_numbers = range(10,100,10)
    rounds=10
    fig = plt.figure()
    ax = fig.add_subplot(111)
    first=True
    for nsteps in step_numbers:
        rms = []
        rms_fluc = []
        std_err = []
        for i in range(rounds):
            repwalk=Reptation_walker(nsteps=nsteps)
            etedist_list, length_list = repwalk.get_multiple_end_to_end_distances(nsteps=nsteps,nwalks=nwalks,avoid=avoid,limited=limited)
            #RMS end-to-end distance
            rms.append(np.sqrt(np.mean(np.square(etedist_list))))
            #RMS fluctuation estimate
            rms_fluc.append(np.sqrt((np.mean(np.square(etedist_list))-np.mean(etedist_list)**2)*nwalks/(nwalks-1)))
            #Standard error estimate
            std_err.append(np.sqrt((np.mean(np.square(etedist_list))-np.mean(etedist_list)**2)/(nwalks-1)))
        chain_lengths = [nsteps*repwalk.r for i in range(rounds)]

        ax.plot(chain_lengths,rms,color="blue")
        ax.plot(chain_lengths,rms_fluc,color="orange")
        ax.plot(chain_lengths,std_err,color="green")
        if first == True:
            plt.legend(["RMS End-to-End Distance","RMS Fluctuation Estimate","Standard Error Estimate"])
            first=False
    ax.set_xlabel("Chain length")

    plt.suptitle("End-to-end distance measures vs chain length \n for 100 rounds of " +str(nwalks)+" reptation walks (self-avoiding)")
    plt.show()

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n + 1)


def plot2D(title,xlabel,ylabel,xlists,ylists,labels_list=None):
    """xlists and ylists can be lists of lists if more than one series should be plotted."""
    print("xlists:",xlists)
    print("ylists:",ylists)
    plt.figure()
    plt.title(str(title))
    plt.xlabel(xlabel)#,fontsize=16)
    plt.ylabel(ylabel)#,fontsize=16)
    if type(xlists[0]) is list:
        # if not len(xlists)==len(ylists)==len(labels_list):
        #     raise ValueError    # xlists, ylists and labels_list must be of same length
        for i in range(len(xlists)):
            plt.plot(xlists[i],ylists[i],label=labels_list[i])
    else:
        plt.plot(xlists,ylists)
    plt.legend()
    # plt.show()


def main():
    gridwalk = Grid_walker()
    #I PDF:
    # gridwalk.plot_multiple_end_to_end_distances(nwalks=50,avoid=False,limited=False,forced=False)
    # gridwalk.plot_multiple_end_to_end_distances(nwalks=1000, avoid=True,limited=False,forced=False)
    # gridwalk.plot_multiple_end_to_end_distances(avoid=True,limited=True,forced=False,nwalks=50)
    # gridwalk.plot_multiple_end_to_end_distances(avoid=True, limited=True, forced=True, nwalks=50)
    # gridwalk.plot_success_rate_vs_nsteps(limited=False)
    # gridwalk.plot_success_rate_vs_nsteps(limited=True)
    # gridwalk.normplot(nwalks=100,avoid=True)
    # gridwalk.plot_multiple_end_to_end_distances_holdon(nwalks=1000)


    # print(gridwalk.get_multiple_end_to_end_distances(nwalks=10,avoid=False))
    # gridwalk.walk_without_avoid(nsteps=100,limited=False)
    # gridwalk.walk_with_self_avoid(nsteps=50,limited=True)
    # gridwalk.plot_the_walk(beads=False)
    # print(gridwalk.success_rate(nsteps=10,limited=True))
    # print(gridwalk.success_rate(nsteps=10,limited=False))
    # gridwalk.plot_success_rate_vs_nsteps()
    # gridwalk.hist_quotient_length_etedist(nwalks=1000)
    # gridwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=False,forced=True)
    # gridwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=True,forced=True)
    # gridwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=True,forced=False)
    # gridwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=False)
    # gridwalk.plot_success_rate_vs_nsteps()

    chainwalk = Freely_jointed_chain()
    #In PDF:
    #chainwalk.plot_bead_size_variation(nsteps=15, nwalks=1000, limited=False, forced=False)
    #chainwalk.plot_bead_size_variation(nsteps=15, nwalks=1000, limited=True, forced=False)
    #chainwalk.plot_bead_size_variation(nsteps=15,nwalks=1000,limited=True,forced=True)

    #chainwalk.plot_multiple_end_to_end_distances(nwalks=50)
    #chainwalk.plot_multiple_end_to_end_distances(nwalks=50,avoid=True,limited=True,forced=False)
    #chainwalk.plot_multiple_end_to_end_distances(nwalks=1000,avoid=True,limited=False,forced=False)
    #chainwalk.plot_multiple_end_to_end_distances(nwalks=50,avoid=True,limited=True,forced=True)
    #chainwalk.plot_bead_size_variation(nwalks=50,limited=True,forced=True)
    #chainwalk.plot_bead_size_variation(nwalks=50,limited=True,forced=False)
    #chainwalk.plot_bead_size_variation(limited=True,forced=True,success=True)
    #chainwalk.plot_bead_size_variation(limited=True,forced=False,success=True)
    #chainwalk.plot_bead_size_variation(limited=False,forced=False,success=True)
    #chainwalk.plot_success_rate_vs_nsteps(limited=False)
    #chainwalk.plot_success_rate_vs_nsteps(limited=True)
    #chainwalk.plot_multiple_end_to_end_distances_holdon(nwalks=1000)

    #chainwalk.normplot(nwalks=1000)
    # chainwalk.plot_multiple_end_to_end_distances(nwalks=2,avoid=True)
    # chainwalk.plot_multiple_end_to_end_distances(nwalks=50,avoid=True,limited=True)
    # chainwalk.walk_without_avoid(nsteps=1000,limited=False)
    # chainwalk.plot_the_walk(beads=False)
    # chainwalk.walk_with_self_avoid(nsteps=20,limited=True)
    # chainwalk.plot_the_walk(beads=False)
    # print(chainwalk.success_rate(nsteps=10,limited=True))
    # print(chainwalk.success_rate(nsteps=10,limited=False))
    # # chainwalk.plot_success_rate_vs_nsteps(limited=True)
    # chainwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=False,forced=True)
    # chainwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=True,forced=True)
    # chainwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=True,forced=False)
    # chainwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=False)
    # chainwalk.plot_etedist_parameters_multiple_versions()
    # chainwalk.plot_success_rate_vs_nsteps()

    # dirwalk = Directed_walker()
    # dirwalk.walk_without_avoid(nsteps=1000)
    # dirwalk.plot_the_walk(beads=False)
    # dirwalk.plot_multiple_end_to_end_distances(avoid=False)

    grid_walker_stepl_variations = Grid_walker_stepl_variations()
    # grid_walker_stepl_variations.walk_without_avoid(nsteps=50)
    # grid_walker_stepl_variations.variate_rho = True
    # grid_walker_stepl_variations.walk_with_self_avoid(nsteps=20,limited=True)
    # grid_walker_stepl_variations.plot_the_walk(beads=True)
    # print(grid_walker_stepl_variations.success_rate(nsteps=10,limited=True))
    # print(grid_walker_stepl_variations.success_rate(nsteps=10,limited=False))
    # grid_walker_stepl_variations.plot_success_rate_vs_nsteps()
    # grid_walker_stepl_variations.hist_quotient_length_etedist(nwalks=1000)
    # grid_walker_stepl_variations.plot_bead_size_variation()
    # grid_walker_stepl_variations.plot_bead_size_variation(success=True)
    # grid_walker_stepl_variations.plot_multiple_end_to_end_distances()
    # grid_walker_stepl_variations.plot_etedist_parameters_multiple_versions()
    #grid_walker_stepl_variations.normplot(nwalks=1000)


    chainwalk_stepl_variations = Freely_jointed_chain_stepl_variations(distribution="N") #TODO: Very bad performance
    #In Pdf:
    #chainwalk_stepl_variations.plot_bead_size_variation(limited=True, forced=True, success=True)
    #chainwalk_stepl_variations.plot_bead_size_variation(limited=True, forced=False, success=True)
    #chainwalk_stepl_variations.plot_bead_size_variation(limited=False, forced=False, success=True)

    # chainwalk_stepl_variations.variate_rho = True
    # chainwalk_stepl_variations.walk_with_self_avoid(nsteps=20, limited=True)
    # chainwalk_stepl_variations.plot_the_walk(beads=True)
    # print(chainwalk_stepl_variations.success_rate(nsteps=10,limited=True))
    # print(chainwalk_stepl_variations.success_rate(nsteps=10,limited=False))
    # chainwalk_stepl_variations.hist_quotient_length_etedist(nwalks=1000)
    # chainwalk_stepl_variations.plot_bead_size_variation(limited=True)
    # chainwalk_stepl_variations.plot_bead_size_variation(success=True,limited=False,nsteps=50) #Bad performace due to self.r too short too often (immediate self-intersection)
    # chainwalk_stepl_variations.plot_etedist_parameters_multiple_versions()
    # chainwalk.plot_bead_size_variation(limited=True,forced=True,success=True)
    # chainwalk_stepl_variations.normplot(nwalks=1000)


    reptationwalk = Reptation_walker(nsteps=10)
    #reptationwalk.walk_with_self_avoid()
    #reptationwalk.plot_the_walk(beads=True)
    reptation_plot_multiple_end_to_end_distances(nwalks=10000)

main()
