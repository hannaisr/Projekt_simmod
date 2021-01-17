import random as rnd
import numpy as np
import re
import matplotlib.pyplot as plt # For 3D plot
from mpl_toolkits.mplot3d import Axes3D # For 3D plot
from scipy.stats import norm # histogram fitting
from scipy.stats import normaltest
from scipy import stats as scipy_stats

class Walker():
    """Walked positions are stored in a list"""
    # Initiate list of visited points
    rho0 = 0.3 # size of first bead, the very least 1/2 step length.
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
    last_r = my

    variate_rho=False

    def walk_one_step(self, limited=False):
        """Implemented by child class"""
        pass

    def walk_without_avoid(self,nsteps=100, limited=False):
        """Walk nsteps steps of non-self-avoiding random walk."""
        self.restart()
        for i in range(nsteps):
            self.walk_one_step(limited)

    def walk_with_self_avoid_forced(self,nsteps=100,limited=True,maxfails=np.inf,avoidLastStep=True):
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
                    try_again=self.test_avoid(avoidLastStep=avoidLastStep)
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

    def walk_with_self_avoid(self,nsteps=100,limited=True,maxfails=np.inf,avoidLastStep=True):
        """Walk nsteps steps of self-avoiding random walk. Returns the number of walks that failed."""
        nfails = 0

        self.restart()
        try_again = True
        # Try to assemble a sequence until successful
        while try_again is True:
            try_again = False
            for i in range(nsteps):
                try_again = self.walk_one_step(limited)
                if try_again is False:
                    # Test if the site is already occupied.
                    try_again=self.test_avoid(avoidLastStep=avoidLastStep)
                # In case of self interception, break attempt immediately
                if try_again is True:
                    # print('Managed to walk',len(self.visited_points)-2,'steps')
                    nfails += 1
                    if nfails >= maxfails:
                        # print("Maximum fails reached")
                        return nfails
                    self.restart()
                    break
        # print('Managed to walk', len(self.visited_points) - 1, 'steps')
        return nfails

    def generate_rho(self):
        if self.variate_rho is True:
            return(self.r-self.visited_points[-1][3]-0.0000001) # Current and last bead radii always cover the whole length of the step
        else:
            return(self.rho0)

    def test_avoid(self,avoidLastStep=True):
        """Implemented by child class"""
        pass

    def restart(self):
        """Resets list of visited points."""
        self.origin = self.origin[:3]
        self.origin.append(self.rho0)
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
        return x,y,z,rho

    def get_success_rate(self,nsteps=20,limited=True,avoidLastStep=True):
        """Calculates success rate from nsuccessful_walks number of successful walks."""
        nfails = 0
        nsuccess = 0
        nwalks = 1000
        while (nfails+nsuccess) < nwalks:
            maxfails = nwalks-(nsuccess+nfails) # Maximum number of failed attempts before breaking loop of self avoiding walk
            nfails += self.walk_with_self_avoid(nsteps=nsteps,limited=limited,maxfails=maxfails,avoidLastStep=avoidLastStep)
            if nfails+nsuccess < nwalks:
                nsuccess += 1
        return nsuccess/(nsuccess+nfails)

    def get_end_to_end_distance(self):
        """Calculate end-to-end distance of already walked walk."""
        last_point = self.visited_points[-1]
        return np.sqrt((last_point[0]-self.origin[0])**2+(last_point[1]-self.origin[1])**2+(last_point[2]-self.origin[2])**2)

    def get_multiple_end_to_end_distances(self,nsteps=100,nwalks=10,avoid=False,limited=True,forced=False,avoidLastStep=True):
        """Returns a list of end-to-end distances and chain lengths for nwalks number of walks of length nsteps"""
        etedist_list = np.zeros(nwalks)
        length_list = np.zeros(nwalks)
        for i in range(nwalks):
            self.restart()
            if avoid is True:
                if forced is True:
                    self.walk_with_self_avoid_forced(nsteps=nsteps,limited=limited,avoidLastStep=avoidLastStep)
                else:
                    self.walk_with_self_avoid(nsteps=nsteps,limited=limited,avoidLastStep=avoidLastStep)
            else:
                self.walk_without_avoid(nsteps=nsteps)
            etedist_list[i] = self.get_end_to_end_distance()
            length_list[i] = self.length
            # print("Finished",i+1,"walks")
        return etedist_list, length_list

    def plot_multiple_end_to_end_distances(self,nwalks=10,nsteps_list=range(5,14,1),avoid=False,limited=False,forced=False,holdon=False,avoidLastStep=True):
        """Plots end-to-end distance RMS, RMS fluctuation and standard error estimate for nwalks walks by chain length"""
        rms=[]
        rms_fluc = []
        std_err = []
        chain_lengths = []
        for nsteps in nsteps_list:
            print(nsteps)
            etedist_list, length_list = self.get_multiple_end_to_end_distances(nsteps=nsteps,nwalks=nwalks,avoid=avoid,limited=limited,forced=forced,avoidLastStep=avoidLastStep)
            #RMS end-to-end distance
            rms.append(np.sqrt(np.mean(np.square(etedist_list))))
            #RMS fluctuation estimate
            rms_fluc.append(np.sqrt((np.mean(np.square(etedist_list))-np.mean(etedist_list)**2)*nwalks/(nwalks-1)))
            #Standard error estimate
            std_err.append(np.sqrt((np.mean(np.square(etedist_list))-np.mean(etedist_list)**2)/(nwalks-1)))
        chain_lengths = [i * self.my for i in nsteps_list]
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

    def plot_multiple_end_to_end_distances_holdon(self,nwalks=100,show=True):
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
        if show==True:
            plt.show()
        return

    def hist_quotient_length_etedist(self,nsteps=100,nwalks=10,avoid=False,limited=True,forced=False,avoidLastStep=True):
        """Plots the quotient between total chain length and end-to-end distance. Returns list of the quotients"""
        etedist_list, length_list = self.get_multiple_end_to_end_distances(nsteps=nsteps,nwalks=nwalks,avoid=avoid,limited=limited,forced=forced,avoidLastStep=avoidLastStep)
        quotients = [etedist_list[i]/length_list[i] for i in range(len(length_list))]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(x=quotients,bins=15)
        ax.vlines(ymin=0,ymax=nwalks/5,x=np.mean(quotients),color='red',label="Mean")
        ax.set_xlabel("End-to-End distance/length")
        plt.legend()
        plt.show()
        return quotients

    def plot_etedist_normal_parameters(self,nwalks=1000,avoid=False,limited=False,forced=False,show=True,avoidLastStep=True):
        """Plots fitted mu and sigma for end-to-end distance vs length of chain. Returns lists of mus, sigmas, lengths and p-values"""
        etedist_lists = []
        nsteps_list = np.arange(2,25,1)
        length_list = list(nsteps_list*self.my)
        mu_list = []
        sigma_list = []
        pval_list = []
        # Get end-to-end distances
        for nsteps in nsteps_list:
            etedist_list, _ = self.get_multiple_end_to_end_distances(nsteps=nsteps,nwalks=nwalks,avoid=avoid,limited=limited,forced=forced,avoidLastStep=avoidLastStep)
            [mu,sigma] = norm.fit(etedist_list)
            mu_list.append(mu)
            sigma_list.append(sigma)
            p=normaltest(etedist_list)
            pval_list.append(p[1])
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
        plot2D(title,"Chain length",None,[length_list,length_list],[mu_list,sigma_list],labels,show=show)
        print("List of mus:",mu_list,"\n List of sigmas:",sigma_list)
        return mu_list, sigma_list, length_list, pval_list

    def plot_etedist_normal_parameters_multiple_methods(self):
        """Plots end-to-end distance mu and sigma for different methods of random walk"""
        mu1,sigma1,lenlist1,pval1 = self.plot_etedist_normal_parameters(limited=False,show=False)
        mu3,sigma3,lenlist3,pval3 = self.plot_etedist_normal_parameters(avoid=True,limited=True,show=False)
        mu4,sigma4,lenlist4,pval4 = self.plot_etedist_normal_parameters(limited=False,forced=True,avoid=True,show=False)
        mu2,sigma2,lenlist2,pval2 = self.plot_etedist_normal_parameters(limited=False,avoid=True,forced=False,show=False)
        plot2D("End-to-end distance mu vs chain length\n for "+str(self.name),"Chain length", "End-to-end distance mu", [lenlist1,lenlist2,lenlist3,lenlist4],[mu1,mu2,mu3,mu4],["Not limited, not avoiding","Self avoiding", "Self avoiding limited","Forced self avoiding"],show=False,scale='log')
        plot2D("End-to-end distance sigma vs chain length\n for "+str(self.name),"Chain length", "End-to-end distance sigma", [lenlist1,lenlist2,lenlist3,lenlist4],[sigma1,sigma2,sigma3,sigma4],["Not limited, not avoiding","Self avoiding", "Self avoiding limited","Forced self avoiding"],show=False,scale='log')
        plot2D("Normal distribution p-value \n for "+str(self.name),"Chain length", "P-value", [lenlist1,lenlist2,lenlist3,lenlist4],[pval1,pval2,pval3,pval4],["Not limited, not avoiding","Self avoiding", "Self avoiding limited","Forced self avoiding"],show=False)
        plt.show()

    def get_spec_distance(self,spec='max'):
        """Calculate maximum distance between points in walked walk
        spec can be 'max', 'median' or 'mean'"""
        dists = []
        counter = 0
        for i in self.visited_points:
            counter += 1
            for j in self.visited_points[counter:]:
                dists.append(np.sqrt((i[0]-j[0])**2+(i[1]-j[1])**2+(i[2]-j[2])**2))
        if spec == 'max':
            return max(dists)
        if spec == 'mean':
            return np.mean(dists)
        if spec == 'median':
            return np.median(dists)

    def get_mean_maximum_distance(self,avoid=False,limited=False,nsteps=15,nwalks=1000,whatiwant="mean",spec="max",avoidLastStep=True):
        """Calculate mean maximum difference for specific walk.
        Possible options for whatiwant:
        'mean' (returns mean value)
        'stdev' (returns standard deviation)"""
        maxdist_list = []
        for i in range(nwalks):
            if avoid is True:
                self.walk_with_self_avoid(limited=limited,nsteps=nsteps,avoidLastStep=avoidLastStep)
            else:
                self.walk_without_avoid(limited=limited,nsteps=nsteps)
            maxdist_list.append(self.get_spec_distance(spec=spec))
        if whatiwant=="mean":
            return np.mean(maxdist_list)
        if whatiwant=="stdev":
            return stdev(maxdist_list)

    def plot_distance_between(self,maxSteps=15,avoid=False,limited=False,nwalks=1000,whatiwant="mean",spec="max"):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ext = ""
        if limited == True:
            ext += ", limited"
        rhos = list(range(1, 5))
        #rhos.append(4.99)#Bad efficiency
        cmap=get_cmap(maxSteps)

        for nsteps in range(2,maxSteps+1):
            print("nsteps=",nsteps)
            qs=[] #Quotients bead size/expected value of step length
            dist=[]
            for rho in rhos:
                print("rho=",rho)
                self.rho0 = rho / 10
                qs.append(self.rho0 / self.my)
                dist.append(self.get_mean_maximum_distance(avoid=avoid,limited=limited,nsteps=nsteps,nwalks=nwalks,whatiwant=whatiwant,spec=spec))
            ax.plot(qs,dist,label="steps= "+str(nsteps),color=cmap(nsteps))
        plt.legend()
        plt.show()

    def plot_success_rate_vs_nsteps(self,step_numbers=range(2,25,2),limited=True,avoidLastStep=True):
        """Gets success rate to number of steps in walk."""
        # step_numbers = range(10,120,10)   # Grid
        # step_numbers = range(1,25,5)        # Freely jointed
        success_rates = []
        for nsteps in step_numbers:
            success_rates.append(self.get_success_rate(nsteps,limited,avoidLastStep=avoidLastStep))
            print(nsteps)
        title = "Success rate vs number of steps for \n"+str(self.name)+", "
        if limited is True:
            title += 'limited'
        else:
            title += 'not limited'
        plot2D(title,'Number of steps','Success rate',step_numbers,success_rates)
        plt.show()
        return step_numbers,success_rates

    def plot_bead_size_variation(self,nsteps=100,nwalks=10,my=True,rhos=np.arange(0.099,0.5,0.025),success=False,limited=False,forced=False,holdon=False,avoidLastStep=True):
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
        # rhos = list(range(1,5))
        # rhos.append(4.99)
        for rho in rhos:
            print(rho)
            self.rho0=rho
            if my==True: #Bead size by expected value
                qs.append(self.rho0/self.my)
            else: #Bead size by standard deviation
                qs.append(self.rho0/self.sigma)
            if success==False:#y=RMS
                etedist_list, length_list = self.get_multiple_end_to_end_distances(nsteps=nsteps, nwalks=nwalks, avoid=True, limited=limited,forced=forced,avoidLastStep=avoidLastStep)
                # RMS end-to-end distance
                rms.append(np.sqrt(np.mean(np.square(etedist_list))))
                # RMS fluctuation estimate
                rms_fluc.append(np.sqrt((np.mean(np.square(etedist_list)) - np.mean(etedist_list) ** 2) * nwalks / (nwalks - 1)))
                # Standard error estimate
                std_err.append(np.sqrt((np.mean(np.square(etedist_list)) - np.mean(etedist_list) ** 2) / (nwalks - 1)))
            else:#y=Success rate
                success_rates.append(self.success_rate(nsteps=nsteps,limited=limited))
        if holdon==True:
            return qs, rms, rms_fluc, std_err, success_rates

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

    def plot_multiple_bead_size_variations_holdon(self,nwalks=1000,nsteps=10,types_of_walks=[1,2,3],data=[0,1,2],rhos=np.arange(0.099,0.5,0.025),plot_rms_grid=False,save=False,avoidLastStep=True,title=None,labels=False,boxes=True,ylabel=None):
        """Plots the relationship between quotient bead size/my or sigma of step length and RMS End-to-End distance or success rate in self-avoiding chains. All chain-types in same plot."""
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel("Bead radius/Expected value of step length")
        ax.set_ylabel(ylabel)
        if title==None:
            title = "End-to-end distance measures vs bead radius \n for " + str(nwalks) + " walks of " + str(nsteps) + "-steps " + self.name
        plt.suptitle(title)
        cmap = get_cmap(n=len(types_of_walks)*len(nsteps))
        if type(nsteps)==int or type(nsteps)==float:
            nsteps = [nsteps]

        j=0
        for n in nsteps:
            # Self avoid
            if 1 in types_of_walks:
                i = j
                qs, rms, rms_fluc, std_err, success_rates = self.plot_bead_size_variation(nwalks=nwalks,nsteps=n,limited=False,rhos=rhos,forced=False,holdon=True, avoidLastStep=avoidLastStep)
                if 0 in data:
                    ax.plot(qs, rms, label="Self avoiding" if labels==False else labels[0], color=cmap(i))
                if 1 in data:
                    ax.plot(qs, rms_fluc, color=cmap(i))
                if 2 in data:
                    ax.plot(qs, std_err, color=cmap(i))

            # Self avoid:limited
            if 2 in types_of_walks:
                i = j+1
                qs, rms, rms_fluc, std_err, success_rates = self.plot_bead_size_variation(nwalks=nwalks,nsteps=n,limited=True,rhos=rhos,forced=False,holdon=True, avoidLastStep=avoidLastStep)
                if 0 in data:
                    ax.plot(qs, rms, label="Self avoiding: limited" if labels==False else labels[1], color=cmap(i))
                if 1 in data:
                    ax.plot(qs, rms_fluc, color=cmap(i))
                if 2 in data:
                    ax.plot(qs, std_err, color=cmap(i))

            # Self avoid:forced
            if 3 in types_of_walks:
                i = j+2
                qs, rms, rms_fluc, std_err, success_rates = self.plot_bead_size_variation(nwalks=nwalks,nsteps=n,limited=True,rhos=rhos,forced=True,holdon=True, avoidLastStep=avoidLastStep)
                if 0 in data:
                    ax.plot(qs, rms, label="Self avoiding: forced" if labels==False else labels[2], color=cmap(i))
                if 1 in data:
                    ax.plot(qs, rms_fluc, color=cmap(i))
                if 2 in data:
                    ax.plot(qs, std_err, color=cmap(i))

            if plot_rms_grid is True:
                # Rms for self avoiding grid walk of different number of steps (0-19) (obtained by walking 10 000 walks for each number of steps)
                rms_grid = [0.0, 1.0, 1.549128787415688, 1.97453792062852, 2.3518928547023563,
                2.6816039976103854, 3.0023324266310016, 3.301999394306425, 3.5805865441293276,
                3.844866707702622, 4.119465984809196, 4.3475510347780855, 4.590250537824706,
                4.815952657574615, 5.017868073196026, 5.224232000973923, 5.449091667424947,
                5.652397721321457, 5.850743542490989, 6.029278563808443]

                ax.plot(qs, [rms_grid[n]]*len(qs), label=str(n)+" steps",color=cmap(j+types_of_walks[0]-1))

            j+=1


        if boxes == True:
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            if 0 in data:
                textstr = "RMS"
                # place a text box in upper left in axes coords
                ax.text(0.90, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                        verticalalignment='top', bbox=props)
            if 1 in data:
                textstr = "RMS\nfluctuation"
                # place a text box in upper left in axes coords
                ax.text(0.90, 0.3, textstr, transform=ax.transAxes, fontsize=10,
                        verticalalignment='top', bbox=props)
            if 2 in data:
                textstr = "Standard error\nestimate"
                # place a text box in upper left in axes coords
                ax.text(0.90, 0.1, textstr, transform=ax.transAxes, fontsize=10,
                        verticalalignment='top', bbox=props)

        print(qs,rms)

        plt.legend()
        if save==True:
            plt.savefig(title,bbox_inches='tight')
        plt.show()
        return

    def normplot(self,nsteps=100,nwalks=10,avoid=False,limited=True,forced=False,avoidLastStep=True):
        """Tests adherence to normal distribution"""
        etedistlist,l = self.get_multiple_end_to_end_distances(nsteps=nsteps,nwalks=nwalks,avoid=avoid,limited=limited,forced=forced,avoidLastStep=avoidLastStep)

        # Calculate quantiles and least-square-fit curve
        (quantiles, values), (slope, intercept, r) = scipy_stats.probplot(etedistlist, dist='norm')

        # plot results
        plt.plot(values, quantiles, 'ob')
        plt.plot(quantiles * slope + intercept, quantiles, 'r')
        plt.xlabel("End-to-end distance")
        plt.ylabel("Quantiles")

        # define ticks
        ticks_perc = [1, 5, 10, 20, 50, 80, 90, 95, 99]

        # transfrom them from precentile to cumulative density
        ticks_quan = [scipy_stats.norm.ppf(i / 100.) for i in ticks_perc]

        # assign new ticks
        plt.yticks(ticks_quan, ticks_perc)

        # show plot
        ext = ", r=" +str(self.rho0)
        if avoid==True:
            ext+=", self-avoiding"
            if limited == True:
                ext += ", limited"
            if forced == True:
                ext += ", forced"

        plt.suptitle("Normal distribution-fit of End-to-end distance \n for " + str(nwalks) + " walks of " +str(nsteps) +" steps for "+ self.name + ext)
        plt.grid()
        plt.show()

class Grid_walker(Walker):
    def __init__(self):
        self.name = "Grid walker"
        self.shortname = "grid"

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
        return False

    def test_avoid(self,avoidLastStep=True):
        """Test if latest site is already occupied. Return True if so, False if not."""
        # if self.r < 2*self.rho:
        #     return True
        if avoidLastStep is True:
            last = -1
        elif avoidLastStep is False:
            last = -2

        if self.visited_points[0][:3] == self.visited_points[-1][:3]:
            return True
        elif any(t == self.visited_points[-1][:3] for t in [i[:3] for i in self. visited_points[:last]]):
            return True
        return False

class Freely_jointed_chain(Walker):
    def __init__(self):
        self.name = 'Freely jointed chain'
        self.shortname = "FJ"

    def walk_one_step(self, limited=False):
        current_pos = self.visited_points[-1][:]
        # Get bead size
        # rho = self.generate_rho()
        rho = self.rho0

        ###---Old version---###
        # # Get a walking direction to star with
        # theta = rnd.uniform(0,np.pi)
        # phi = rnd.uniform(0,2*np.pi)
        #
        # if limited == True and self.last_direction != 0:
        #     # Define direction to walk back the same way
        #     theta_back = np.pi-self.last_direction[0]
        #     phi_back = np.pi+self.last_direction[1]
        #
        #     while (self.r*np.sin(theta)*np.cos(phi)-self.r*np.sin(theta_back)*np.cos(phi_back))**2+(self.r*np.sin(theta)*np.sin(phi)-self.r*np.sin(theta_back)*np.sin(phi_back))**2+(self.r*np.cos(theta)-self.r*np.cos(theta_back))**2 < (current_pos[3]+rho)**2:
        #        theta = rnd.uniform(0,np.pi)
        #        phi = rnd.uniform(0,2*np.pi)
        #
        # self.last_direction = [theta,phi]
        # # Update the coordinates
        # current_pos[0] += self.r*np.sin(theta)*np.cos(phi)  # x
        # current_pos[1] += self.r*np.sin(theta)*np.sin(phi)  # y
        # current_pos[2] += self.r*np.cos(theta)              # z
        # current_pos[3] = rho
        ###---End of old version---###

        ###---Suggestion---###
        if limited == True and self.last_direction != 0:
            # Define direction to walk back the same way
            theta_back = np.pi-self.last_direction[0]
            phi_back = np.pi+self.last_direction[1]

            #In local coordinate system with "going back" vector as z axis:
            #Generate new bead center
            alpha = np.arccos((self.r ** 2 + self.last_r ** 2 - (2*self.visited_points[-2][3]) ** 2) / (2 * self.r * self.last_r))
            z_rnd = rnd.uniform(-1,np.cos(alpha))
            theta = np.arccos(z_rnd)
            phi = rnd.uniform(0, np.pi * 2)

            #Translate bead center position into global coordinate system.
            # First: define the local coordinate system in terms of the global
            z_hat = [np.sin(theta_back) * np.cos(phi_back), np.sin(theta_back) * np.sin(phi_back),np.cos(theta_back)]  # vector between the most recent two beads
            y_hat = [-np.sin(phi_back), np.cos(phi_back), 0] # Normalized and orthogonal to z_hat
            x_hat = [-np.cos(theta_back) * np.cos(phi_back),-np.cos(theta_back) * np.sin(phi_back), np.sin(theta_back)] # Normalized and orthogonal to z_hat and y hat

            #print(sum([y_hat[i] * y_hat[i] for i in range(3)]))
            #Second: project the bead center position onto the global axes and translate it to origo
            current_pos_origo = [self.r*(np.cos(theta)*z_hat[i]+np.sin(theta)*np.cos(phi)*y_hat[i]+np.sin(theta)*np.sin(phi)*x_hat[i]) for i in range(3)]
            current_pos = [current_pos[i] + current_pos_origo[i] for i in range(3)]
            current_pos.append(rho)

            # vector_length = np.sqrt(sum([current_pos[i]**2 for i in range(3)]))
            # Define direction
            self.last_direction = [np.arctan(current_pos_origo[1]/current_pos_origo[0]), np.arccos(current_pos_origo[2]/self.r)]
        else:
            z_rnd = rnd.uniform(-1, 1)
            theta = np.arccos(z_rnd)
            phi = rnd.uniform(0, np.pi * 2)

            self.last_direction = [theta,phi]
            # Update the coordinates
            current_pos[0] += self.r*np.sin(theta)*np.cos(phi)  # x
            current_pos[1] += self.r*np.sin(theta)*np.sin(phi)  # y
            current_pos[2] += self.r*z_rnd                      # z
            current_pos[3] = rho
            ###---End of Suggestion---###

        # Update list of visited points
        self.visited_points.append(current_pos)
        self.length += self.r
        self.last_r = self.r
        return False

    def test_limited_option(self,limited=True):
        """Tests whether the limited-option actually limits options to > 2*rho from the neighbor of the most recent point"""
        def walk_options_help(limited=True):
            options=[]
            current_pos = self.visited_points[-1][:]
            # Get bead size
            rho = self.rho0
            if limited == True and self.last_direction != 0: #The last step

                # Define direction to walk back the same way
                theta_back = np.pi-self.last_direction[0]
                phi_back = np.pi+self.last_direction[1]
                #In local coordinate system with "going back" vector as z axis:
                #Generate new bead center
                self.alpha = np.arccos((self.r ** 2 + self.last_r ** 2 - (2*self.visited_points[-2][3]) ** 2) / (2 * self.r * self.last_r))
                #TESTING STEP: Generate many options for the position of the next bead. All options are stored in
                for i in range(1000):
                    z_rnd = rnd.uniform(-1, np.cos(self.alpha))
                    phi = rnd.uniform(0, np.pi * 2)
                    theta = np.arccos(z_rnd)

                    #Translate bead center position into global coordinate system.
                    # First: define the local coordinate system in terms of the global
                    z_hat = [np.sin(theta_back) * np.cos(phi_back), np.sin(theta_back) * np.sin(phi_back),np.cos(theta_back)]  # vector between the most recent two beads
                    y_hat = [-np.sin(phi_back), np.cos(phi_back), 0] #Orthogonal
                    x_hat = [-np.cos(theta_back) * np.cos(phi_back),-np.cos(theta_back) * np.sin(phi_back), np.sin(theta_back)] #Orthogonal
                    #Second: project the bead center position onto the global axes and translate it to origo
                    option_direction = [self.r*(np.cos(theta)*z_hat[i]+np.sin(theta)*np.cos(phi)*y_hat[i]+np.sin(theta)*np.sin(phi)*x_hat[i]) for i in range(3)]
                    option = [current_pos[i] + option_direction[i] for i in range(3)]
                    options.append(option)
                current_pos=option
                current_pos.append(self.rho0)
            elif limited==False and self.last_direction != 0:#The last step

                for i in range(1000):
                    option=[0,0,0,self.rho0]
                    z_rnd = rnd.uniform(-1, 1)
                    phi = rnd.uniform(0, np.pi * 2)
                    theta = np.arccos(z_rnd)
                    option[0] = current_pos[0]+self.r * np.sin(theta) * np.cos(phi)  # x
                    option[1] = current_pos[1]+self.r * np.sin(theta) * np.sin(phi)  # y
                    option[2] = current_pos[2]+z_rnd  # z
                    options.append(option)
                current_pos = option
                current_pos.append(self.rho0)
            else: #The first step
                z_rnd = rnd.uniform(-1, 1)
                phi = rnd.uniform(0, np.pi * 2)
                theta = np.arccos(z_rnd)

                self.last_direction = [theta,phi]
                # Update the coordinates
                current_pos[0] += self.r*np.sin(theta)*np.cos(phi)  # x
                current_pos[1] += self.r*np.sin(theta)*np.sin(phi)  # y
                current_pos[2] += z_rnd                             # z
                current_pos[3] = rho
            # Update list of visited points
            self.visited_points.append(current_pos)
            self.length += self.r
            self.last_r = self.r
            return options

        #Generate the 3-pearl (2-steps) limited walk
        self.restart()
        for i in range(2):
            options = walk_options_help(limited)
        #Check to see the distance between the center of the first pearl and each option
        dist=[]
        last=self.visited_points[-3]
        for opt in options:
            dist.append(np.sqrt(sum([(opt[i]-last[i])**2 for i in range(3)])))
        print("Distance between pearl 1 and options for pearl 3:\nmin: ",min(dist), ", 2*rho: ",2*self.rho0,", max: ",max(dist),", r:",self.r)

    def test_avoid(self,avoidLastStep=True):
        """Test if latest site is already occupied - continuous case"""
        cp = self.visited_points[-1]    # Current position
        # print(cp)
        #The distance between successive sphere centres is self.r. Interception between any two spheres occurs if their centres are less apart than their diameter
        # if self.r < 2*self.rho:
        #     return True
        if avoidLastStep is True:
            last = -1
        elif avoidLastStep is False:
            last = -2

        for point in self.visited_points[:last]:
            r_centres = np.sqrt((point[0] - cp[0])**2 + (point[1] - cp[1])**2 + (point[2] - cp[2])**2)
            if r_centres < point[3]+cp[3]:
                # Self-intercept - needs to restart the process
                return True
        return False
        # TODO Implement method for avoiding previous _paths_, not just previous positions.

class Reptation_walker(Grid_walker):
    """Reptation algorithm to generate SAWs. Algoritm on hiskp.uni-bonn... pg. 4"""
    origin = [0,0,0]
    def __init__(self,nsteps=100,name='Reptation walker'):
        # Generate a basis self-avoiding walk
        self.gridwalk = Grid_walker()
        self.gridwalk.walk_with_self_avoid(nsteps=nsteps)
        self.visited_points = self.gridwalk.visited_points #The gridwalk.visited_points will be modified as saw is
        self.nsteps=nsteps
        self.name = name
        self.shortname = "rept"

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
        self.shortname = "dirw"

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
        self.name = "Grid walker with step length variation"
        self.shortname = "steplvar grid"

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
        self.name = "Freely jointed chain with step length variations"
        self.shortname = "FJ, " + distribution + " varied stepl"

    def walk_one_step(self, limited=False):
        if self.distribution == "N":
            self.r = rnd.normalvariate(self.my,self.sigma)  # Varation in step length
        elif self.distribution == "exp":
            self.r = rnd.expovariate(1/self.my)  # Varation in step length
        if self.r+self.last_r <= 2*self.rho0:
            # print(self.r,self.last_r)
            self.last_r = self.r
            return True
        return Freely_jointed_chain.walk_one_step(self,limited)

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


def plot2D(title,xlabel,ylabel,xlists,ylists,labels_list=None,scale='linlin',show=True,save=False,fileName=None,labelposition="inside"):
    """xlists and ylists can be lists of lists if more than one series should be plotted."""
    print("xlists:",xlists)
    print("ylists:",ylists)
    plt.figure()
    plt.title(title,fontsize=16)
    plt.xlabel(xlabel,fontsize=16)
    plt.ylabel(ylabel,fontsize=16)
    # Make sure there is enough x-values to plot
    if type(ylists[0])==list and type(xlists[0])!=list:
        n = len(ylists)
        xlists = [xlists]*n
    # Make the plot
    if type(xlists[0]) is list:
        for i in range(len(xlists)):
            plt.plot(xlists[i],ylists[i],label=labels_list[i])
    else:
        plt.plot(xlists,ylists)
    # Generate labels
    if labels_list:
        if labelposition=="inside":
            plt.legend(fontsize=14)
        if labelposition=="outside":
            plt.legend(bbox_to_anchor=(1.05,1),loc="upper left",fontsize=14)
    # Alter scale on axes
    if scale[:3] == 'log':
        plt.xscale('log')
    if scale[-3:] == 'log':
        plt.yscale('log')
    # Change font size for ticks
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # Save file
    if save is True:
        if fileName:
            if fileName[-4:] == ".png":
                name = str(fileName).replace(" ","_")
            else:
                name = str(fileName).replace(" ","_")+".png"
        else:
            name = str(title).replace(" ","_")+".png"
        # print(name)
        plt.savefig(name,bbox_inches='tight')
    # Show plot
    if show is True:
        plt.show()

def plot_success_rate_vs_onXAxis(instances,onXAxis='nsteps',limited=True,bothLimitedAndNot=False,nsteps_list=range(0,25,1),m_list=0.4,show=True,save=False,scale='linlin',labelposition="inside",avoidLastStep=True,title=None,labels_list=None):
    """Plots success rate vs number of steps in walk for various instances, or just one. Also compares limited with not limited if bothLimitedAndNot is True.
    nsteps is the number of steps in each walk
    m = rho/r_mean (m_list is a list of these values to be tested)
    onXAxis = 'nsteps','bead radius', 'bead volume' or 'bead surface area'"""
    instances = list(instances) # require list
    success_rates = []
    if bothLimitedAndNot is True:
        limited = True
    if not labels_list:
        generate_labels = True
        labels_list = []
    else:
        generate_labels = False
    if type(m_list)==float or type(m_list)==int:
        m_list = [m_list]
    if type(scale) == str:
        scale = [scale]

    if onXAxis=='nsteps':
        X = x_list = list(nsteps_list)
        comparisons = m_list
        name_endings = 'number of steps'
        xlabel = "Number of steps in walk"
    else:
        X = list(m_list)
        comparisons = nsteps_list
        name_endings = onXAxis
        if onXAxis == 'bead radius':
            x_list = list(m_list)
        elif onXAxis == 'bead volume':
            x_list = list(m_list**3*4/3*np.pi)
        elif onXAxis == 'bead surface area':
            x_list = list(m_list**2*4*np.pi)
        xlabel = onXAxis.capitalize()


    # Make plot title if no title has been provided
    if not title:
        title = "Success rate vs "+name_endings
        if len(instances)==1:
            title+=",\n"+instances[0].name
        if bothLimitedAndNot is False:
            if limited is True:
                title+=",\n limited walk"
            else:
                title+=",\n regular walk"
        if avoidLastStep is False:
            title += ",\n no collision between consecutive steps"
        if len(comparisons)==1:
            if onXAxis=='nsteps':
                title+='\n '+str(comparisons)+' steps'
            else:
                title+='\n rho/r='+str(comparisons)
                title+=str(comparisons)
        title = re.sub('\[*\]*','',title)
    fileName = "SR vs "+onXAxis


    # Make labels and store success rates
    for instance in instances:
        for c in comparisons:
            if onXAxis=='nsteps':    # Update rho0 if the values of this is what should be compared
                instance.rho0 = c*instance.my
            else:
                nsteps = c

            # --- label generation ---
            if generate_labels is True:
                label = str()
                if len(instances)!=1:
                    label += ", "+instance.shortname
                if bothLimitedAndNot is True:
                    label += ", lim"
                if len(comparisons)>1:
                    if onXAxis=='nsteps':
                        label += ", rho/r="+str(round(c,3))
                    else:
                        label+=", "+str(c)+" steps"
                if len(label) > 2:
                    label = label[2:]   # Remove the comma and space at the beginning
                labels_list.append(label)
            # --- end of label generation ---

            # Get success rates
            SR = []
            for x in X:
                if onXAxis=='nsteps':
                    nsteps=x
                else:
                    instance.rho0=x*instance.my
                SR.append(instance.get_success_rate(nsteps=nsteps,limited=limited,avoidLastStep=avoidLastStep))
                # print(nsteps)
            success_rates.append(SR)

            if bothLimitedAndNot is True:
                # --- Label generation ---
                label = str()
                if len(instances)>1:
                    label += ", "+instance.shortname
                label += ", reg"
                if len(comparisons)>1:
                    if onXAxis=='nsteps':
                        label += ", rho/r="+str(round(c,3))
                    else:
                        label+=", "+str(c)+" steps"
                label = label[2:]
                labels_list.append(label)
                # --- end of label generation ---

                # Success rate
                SR = []
                for x in X:
                    if onXAxis=='nsteps':
                        nsteps=onXAxis
                    else:
                        instance.rho0=x*instance.my
                    SR.append(instance.get_success_rate(nsteps=nsteps,limited=False,avoidLastStep=avoidLastStep))
                    # print(nsteps)
                success_rates.append(SR)
            print(instance.name)

    # Plot
    for instance in instances:
        fileName += " " + instance.shortname
    if limited==True:
        fileName += " lim"
    if limited==False or bothLimitedAndNot==True:
        fileName += " reg"
    fileName +=" m" + str(m_list) + " nsteps"+str(nsteps_list)
    if avoidLastStep is True:
        fileName += " avoidLastStep"
    fileName = re.sub('\W+','_',fileName)
    for sc in scale:
        newFileName=fileName+"_"+sc
        plot2D(title,xlabel,"Success rate",x_list,success_rates,labels_list,scale=sc,show=show,fileName=newFileName,save=save,labelposition=labelposition)
    return x_list, success_rates

def plot_success_rate_vs_r(instances,nsteps=10,limited=True,bothLimitedAndNot=False,r_range=np.arange(1,10,1),m_list=0.4,show=True,save=False,scale='linlin',labelposition="inside",avoidLastStep=True):
    """Plots success rate vs bead size for one or more instances. Also compares limited with not limited if bothLimitedAndNot is True.
    m_list is bead size/step length, can be a list"""
    instances = list(instances) # require list
    success_rates = []
    if type(m_list) == float:
        m_list = [m_list]
    if bothLimitedAndNot is True:
        limited = True

    # Make plot title
    title = "Success rate vs step length, " + str(nsteps) + " steps, "
    if avoidLastStep is True:
        title += "\n avoiding last step"
    else:
        title += "\n not avoiding last step"

    # Make labels and store success rates
    labels_list = []
    for instance in instances:
        for m in m_list:
            #  Labels
            if (limited is True):
                labels_list.append(str(instance.shortname)+", limited, rho/r="+str(m))
            else:
                labels_list.append(str(instance.shortname)+", rho/r="+str(m))
            # Success rates
            SR = []
            for r in r_range:
                instance.my = r
                instance.rho0 = r*m
                SR.append(instance.get_success_rate(nsteps,limited,avoidLastStep=avoidLastStep))
                print(r)
            success_rates.append(SR)
            if bothLimitedAndNot is True:
                labels_list.append(str(instance.shortname)+", regular, rho/r="+str(m))
                SR = []
                for r in r_range:
                    instance.my = r
                    instance.rho0 = r*m
                    SR.append(instance.get_success_rate(nsteps,limited=False,avoidLastStep=avoidLastStep))
                    print(r)
                success_rates.append(SR)
            print(instance.name)

    # Plot
    fileName = "SR_vs_r_m" + str(m_list)
    for instance in instances:
        fileName += " " + instance.shortname
    if bothLimitedAndNot is True:
        fileName += " lim and reg "
    elif limited is True:
        fileName += " limited "
    else:
        fileName += " regular "
    print(fileName)
    fileName = re.sub('\W+',' ',fileName) + " " + str(scale)
    if avoidLastStep is True:
        fileName += " avoidLastStep"
    print(fileName)
    plot2D(title,"Step length", "Success rate", list(r_range), success_rates, labels_list,scale=scale,show=show,fileName=fileName,save=save,labelposition=labelposition)
    return r_range, success_rates

def plot_multiple_end_to_end_distances(instances,types_of_walks=[0,1,2,3],data=[0,1,2],nwalks=1000,nsteps_list=range(5,14,1),avoidLastStep=True,labelposition="inside",show=True,save=False):
    """types_of_walks:
            0 = non-self avoiding (regular)
            1 = self-avoiding
            2 = self-avoiding, limited
            3 = self-avoiding, forced
        data:
            0 = rms
            1 = rms_fluc
            2 = std_err"""
    instances = list(instances)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("Chain length")
    title="End-to-end distance measures vs chain length, \n" +str(nwalks)+" walks"
    if len(instances)==1 and instances[0].shortname!='grid':
        title += ', bead radius '+ str(instances[0].rho0)
    else:
        title = title[:-2]
    plt.suptitle(title)
    cmap=get_cmap(n=len(types_of_walks)*len(instances))

    j = 0
    for instance in instances:
        #Non-self avoid
        if 0 in types_of_walks:
            i=j
            chain_lengths, rms, rms_fluc, std_err = instance.plot_multiple_end_to_end_distances(nsteps_list=nsteps_list,nwalks=nwalks,avoid=False,limited=False,forced=False,holdon=True,avoidLastStep=avoidLastStep)
            if 0 in data:
                ax.plot(chain_lengths, rms, label="Non self avoiding, "+instance.shortname,color=cmap(i))
            if 1 in data:
                ax.plot(chain_lengths, rms_fluc,label="Non self avoiding, "+instance.shortname if not 0 in data else None, color=cmap(i))
            if 2 in data:
                ax.plot(chain_lengths, std_err,label="Non self avoiding, "+instance.shortname if not 0 in data and not 1 in data else None,color=cmap(i))

        #Self avoid
        if 1 in types_of_walks:
            i=j+1
            chain_lengths, rms, rms_fluc, std_err = instance.plot_multiple_end_to_end_distances(nsteps_list=nsteps_list,nwalks=nwalks,avoid=True,limited=False,forced=False,holdon=True,avoidLastStep=avoidLastStep)
            if 0 in data:
                ax.plot(chain_lengths, rms, label="Self avoiding, "+instance.shortname,color=cmap(i))
            if 1 in data:
                ax.plot(chain_lengths, rms_fluc,label="Self avoiding, "+instance.shortname if not 0 in data else None,color=cmap(i))
            if 2 in data:
                ax.plot(chain_lengths, std_err,label="Self avoiding, "+instance.shortname if not 0 in data and not 1 in data else None, color=cmap(i))
            print(rms)

        #Self avoid:limited
        if 2 in types_of_walks:
            i=j+2
            chain_lengths, rms, rms_fluc, std_err = instance.plot_multiple_end_to_end_distances(nsteps_list=nsteps_list,nwalks=nwalks,avoid=True,limited=True,forced=False,holdon=True,avoidLastStep=avoidLastStep)
            if 0 in data:
                ax.plot(chain_lengths, rms, label="Self avoiding: limited, "+instance.shortname,color=cmap(i))
            if 1 in data:
                ax.plot(chain_lengths, rms_fluc,label="Self avoiding: limited, "+instance.shortname if not 0 in data else None, color=cmap(i))
            if 2 in data:
                ax.plot(chain_lengths, std_err,label="Self avoiding: limited, "+instance.shortname if not 0 in data and not 1 in data else None,color=cmap(i))

        #Self avoid:forced
        if 3 in types_of_walks:
            i=j+3
            chain_lengths, rms, rms_fluc, std_err = instance.plot_multiple_end_to_end_distances(nsteps_list=nsteps_list,nwalks=nwalks,avoid=True,limited=True,forced=True,holdon=True,avoidLastStep=avoidLastStep)
            if 0 in data:
                ax.plot(chain_lengths, rms, label="Self avoiding: forced, "+instance.shortname,color=cmap(i))
            if 1 in data:
                ax.plot(chain_lengths, rms_fluc,label="Self avoiding: forced, "+instance.shortname if not 0 in data else None,color=cmap(i))
            if 2 in data:
                ax.plot(chain_lengths, std_err,label="Self avoiding: forced, "+instance.shortname if not 0 in data and not 1 in data else None, color=cmap(i))
        j+=1

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    if 0 in data:
        textstr="RMS"
        # place a text box in upper left in axes coords
        ax.text(0.90, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
    if 1 in data:
        textstr = "RMS\nfluctuation"
        # place a text box in upper left in axes coords
        ax.text(0.90, 0.3, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
    if 2 in data:
        textstr = "Standard error\nestimate"
        # place a text box in upper left in axes coords
        ax.text(0.90, 0.1, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)

    if labelposition=="inside":
        plt.legend()
    elif labelposition=="outside":
        plt.legend(bbox_to_anchor=(1.05,1),loc="upper left")
    if save==True:
        plt.savefig(title+str(types_of_walks)+str(data),bbox_inches='tight')
    if show==True:
        plt.show()
    return


def main():
    gridwalk = Grid_walker()
    #-- I PDF:
    #-RMS measures:
    # gridwalk.plot_multiple_end_to_end_distances(nwalks=50,avoid=False,limited=False,forced=False)
    # gridwalk.plot_multiple_end_to_end_distances(nwalks=1000, avoid=True,limited=False,forced=False)
    # gridwalk.plot_multiple_end_to_end_distances(avoid=True,limited=True,forced=False,nwalks=50)
    # gridwalk.plot_multiple_end_to_end_distances(avoid=True, limited=True, forced=True, nwalks=50)
    #-RMS measures: self avoiding options on same plot
    # gridwalk.plot_success_rate_vs_nsteps(limited=False)
    #-Normplot
    # gridwalk.normplot(nwalks=100,avoid=True)

    #--Not in PDF
    #-Single walk:
    # gridwalk.walk_without_avoid(nsteps=100,limited=False)
    # gridwalk.walk_with_self_avoid(nsteps=50,limited=True)
    # gridwalk.plot_the_walk(beads=False)
    # print(gridwalk.get_success_rate(nsteps=10,limited=True))
    # print(gridwalk.get_success_rate(nsteps=10,limited=False))
    #-ETE distance histogram
    # gridwalk.hist_quotient_length_etedist(nwalks=1000)
    # gridwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=False,forced=True)
    # gridwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=True,forced=True)
    # gridwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=True,forced=False)
    # gridwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=False)
    # gridwalk.plot_success_rate_vs_nsteps()

    chainwalk = Freely_jointed_chain()
    #--In PDF:
    #-By bead size
    #chainwalk.plot_bead_size_variation(nsteps=15, nwalks=1000, limited=False, forced=False)
    #chainwalk.plot_bead_size_variation(nsteps=15, nwalks=1000, limited=True, forced=False)
    #chainwalk.plot_bead_size_variation(nsteps=15,nwalks=1000,limited=True,forced=True)
    #chainwalk.plot_multiple_end_to_end_distances_holdon(nwalks=1000)

    #-RMS measures
    #chainwalk.plot_multiple_end_to_end_distances(nwalks=50)
    #chainwalk.plot_multiple_end_to_end_distances(nwalks=50,avoid=True,limited=True,forced=False)
    #chainwalk.plot_multiple_end_to_end_distances(nwalks=1000,avoid=True,limited=False,forced=False)
    #chainwalk.plot_multiple_end_to_end_distances(nwalks=50,avoid=True,limited=True,forced=True)
    #-Success rate
    #chainwalk.plot_success_rate_vs_nsteps(limited=False)
    #chainwalk.plot_success_rate_vs_nsteps(limited=True)
    #-Single walk
    # chainwalk.walk_without_avoid(nsteps=1000,limited=False)
    # chainwalk.plot_the_walk(beads=False)
    # chainwalk.walk_with_self_avoid(nsteps=20,limited=True,avoidLastStep=False)
    # chainwalk.plot_the_walk(beads=True)
    # chainwalk.test_limited_option(limited=True)

    #-Other
    # print(chainwalk.get_success_rate(nsteps=10,limited=True))
    # print(chainwalk.get_success_rate(nsteps=10,limited=False))
    # chainwalk.plot_success_rate_vs_nsteps(limited=True)
    # chainwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=False)
    # chainwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=False,forced=True)
    # chainwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=True,forced=True)
    # chainwalk.hist_quotient_length_etedist(nsteps=15,nwalks=1000,avoid=True,limited=True,forced=False)
    # chainwalk.plot_etedist_parameters_multiple_versions()

    # dirwalk = Directed_walker()
    # dirwalk.walk_without_avoid(nsteps=1000)
    # dirwalk.plot_the_walk(beads=False)
    # dirwalk.plot_multiple_end_to_end_distances(avoid=False)

    grid_walker_stepl_variations = Grid_walker_stepl_variations()
    # grid_walker_stepl_variations.walk_without_avoid(nsteps=50)
    # grid_walker_stepl_variations.variate_rho = True
    # grid_walker_stepl_variations.walk_with_self_avoid(nsteps=20,limited=True)
    # grid_walker_stepl_variations.plot_the_walk(beads=True)
    # print(grid_walker_stepl_variations.get_success_rate(nsteps=10,limited=True))
    # print(grid_walker_stepl_variations.get_success_rate(nsteps=10,limited=False))
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
    # chainwalk_stepl_variations.walk_with_self_avoid(nsteps=15, limited=True)
    # chainwalk_stepl_variations.plot_the_walk(beads=True)
    # print(chainwalk_stepl_variations.get_success_rate(nsteps=10,limited=True))
    # print(chainwalk_stepl_variations.get_success_rate(nsteps=10,limited=False))
    # chainwalk_stepl_variations.hist_quotient_length_etedist(nwalks=1000)
    # chainwalk_stepl_variations.plot_bead_size_variation(limited=True)
    # chainwalk_stepl_variations.plot_bead_size_variation(success=True,limited=False,nsteps=50) #Bad performace due to self.r too short too often (immediate self-intersection)
    # chainwalk_stepl_variations.plot_etedist_parameters_multiple_versions()
    # chainwalk.plot_bead_size_variation(limited=True,forced=True,success=True)
    # chainwalk_stepl_variations.normplot(nwalks=1000)


    reptationwalk = Reptation_walker(nsteps=10)
    #reptationwalk.walk_with_self_avoid()
    #reptationwalk.plot_the_walk(beads=True)
    # reptation_plot_multiple_end_to_end_distances(nwalks=10000)


    ### SUCCESS RATE PLOTS ###

    scales = ['linlin','linlog']

    # ==========================================================================
    #-------------------------------------#
    # # SUCCESS RATE VS NUMBER OF STEPS # #
    #-------------------------------------#
    nsteps_list = range(1,16,1)
    ## FREELY JOINTED, AVOID LAST BEAD
    m_list1 = np.arange(0.099,0.51,0.1)
    # Limited, lin-lin and lin-log scale, comparison multiple rho/r
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='nsteps',limited=True,nsteps_list=nsteps_list,m_list=m_list1,show=False,save=True,scale=scales)
    # Regular, lin-lin and lin-log scale, comparison multiple rho/r
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='nsteps',limited=False,nsteps_list=nsteps_list,m_list=m_list1,show=False,save=True,scale=scales)
    # Complarison limited and regular
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='nsteps',bothLimitedAndNot=True,nsteps_list=nsteps_list,m_list=0.499,show=False,save=True,scale=scales)

    # ## FREELY JOINTED, DON'T AVOID LAST BEAD
    m_list2 = np.arange(0.1,1.0,0.1)
    # Limited, lin-lin and lin-log scale, comparison multiple rho/r
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='nsteps',limited=True,nsteps_list=nsteps_list,m_list=m_list2,show=False,save=True,scale=scales,avoidLastStep=False,labelposition="outside")
    # Regular, lin-lin and lin-log scale, comparison multiple rho/r
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='nsteps',limited=False,nsteps_list=nsteps_list,m_list=m_list2,show=False,save=True,scale=scales,avoidLastStep=False,labelposition="outside")
    # Complarison limited and regular
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='nsteps',bothLimitedAndNot=True,nsteps_list=nsteps_list,m_list=0.499,show=False,save=True,scale=scales,avoidLastStep=False)

    # ## FREELY JOINTED NORMAL VARIATED, AVOID LAST BEAD
    m_list = m_list1
    # Limited, lin-lin and lin-log scale, comparison multiple rho/r
    plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='nsteps',limited=True,nsteps_list=nsteps_list,m_list=m_list,show=False,save=True,scale=scales)  # Limited is a (somewhat) straight line
    # Regular, lin-lin and lin-log scale, comparison multiple rho/r
    plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='nsteps',limited=False,nsteps_list=nsteps_list,m_list=m_list,show=False,save=True,scale=scales)  # Limited is a (somewhat) straight line
    # Complarison limited and regular
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='nsteps',bothLimitedAndNot=True,nsteps_list=nsteps_list,m_list=0.499,show=False,save=True,scale=scales)


    ## FREELY JOINTED NORMAL VARIATED, DON'T AVOID LAST BEAD
    # Limited, lin-lin and lin-log scale, comparison multiple rho/r
    plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='nsteps',limited=True,nsteps_list=nsteps_list,m_list=np.arange(0.1,0.6,0.05),show=False,save=True,scale=scales,avoidLastStep=False,labelposition="outside")
    # Regular, lin-lin and lin-logscale, comparison multiple rho/r
    plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='nsteps',limited=False,nsteps_list=nsteps_list,m_list=np.arange(0.1,0.6,0.05),show=False,save=True,scale=scales,avoidLastStep=False,labelposition="outside")
    # # Complarison limited and regular
    # plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='nsteps',bothLimitedAndNot=True,nsteps_list=nsteps_list,m_list=0.499,show=False,save=True,scale=scales,avoidLastStep=False)

    ## GRID
    # Limited, lin-lin and lin-log scale
    plot_success_rate_vs_onXAxis([gridwalk],onXAxis='nsteps',limited=True,nsteps_list=range(1,100,1),show=False,save=True,scale=scales)
    # Regular, lin-lin and lin-log scale
    plot_success_rate_vs_onXAxis([gridwalk],onXAxis='nsteps',limited=False,nsteps_list=range(1,25,1),show=False,save=True,scale=scales)
    # Comparison between limited and regular, lin-lin and lin-log scale
    plot_success_rate_vs_onXAxis([gridwalk],onXAxis='nsteps',bothLimitedAndNot=True,nsteps_list=range(1,25,1),show=False,save=True,scale=scales)

    ## COMPARISON OF ALL THREE METHODS
    plot_success_rate_vs_onXAxis([gridwalk,chainwalk,chainwalk_stepl_variations],onXAxis='nsteps',bothLimitedAndNot=True,nsteps_list=nsteps_list,m_list=0.499,show=False,save=True,scale=scales,labelposition="outside")

    #---------------------------------#
    # # SUCCESS RATE VS BEAD RADIUS # #
    #---------------------------------#
    # FREELY JOINTED, AVOID LAST BEAD
    m_list1=np.arange(0.099,0.5,0.05)
    # Limited, comparsion of different numbers of steps
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='bead radius',limited=True,nsteps_list=np.arange(2,15,1),m_list=m_list1,show=False,save=True,labelposition="outside",avoidLastStep=True)
    # Regular, comparsion of different numbers of steps
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='bead radius',limited=False,nsteps_list=np.arange(2,15,1),m_list=m_list1,show=False,save=True,labelposition="outside",avoidLastStep=True)
    # Comparison between regular and limited, 10 steps
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='bead radius',bothLimitedAndNot=True,nsteps_list=[10],m_list=m_list1,show=False,save=True,labelposition="outside",avoidLastStep=True,title = "Comparison between success rates for limited and regular walk,\n freely jointed chain")

    # FREELY JOINTED, DON'T AVOID LAST BEAD
    m_list2=np.arange(0.099,1.0,0.05)
    # Limited, comparsion of different numbers of steps
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='bead radius',limited=True,nsteps_list=np.arange(2,15,1),m_list=m_list2,show=False,save=True,labelposition="outside",avoidLastStep=False)
    # Regular, comparsion of different numbers of steps
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='bead radius',limited=False,nsteps_list=np.arange(2,15,1),m_list=m_list2,show=False,save=True,labelposition="outside",avoidLastStep=False)
    # Comparison between regular and limited, 10 steps
    plot_success_rate_vs_onXAxis([chainwalk],onXAxis='bead radius',bothLimitedAndNot=True,nsteps_list=[10],m_list=m_list2,show=False,save=True,labelposition="outside",avoidLastStep=False,title = "Comparison between success rates for limited and regular walk,\n freely jointed chain,\n no collision between consecutive steps")

    # FREELY JOINTED NORMAL VARIATED, AVOID LAST BEAD
    m_list=m_list1 #m_list=np.arange(0.099,0.5,0.00)
    # Limited, comparsion of different numbers of steps
    plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='bead radius',limited=True,nsteps_list=np.arange(2,15,1),m_list=m_list,show=False,save=True,labelposition="outside",avoidLastStep=True)
    # Regular, comparsion of different numbers of steps
    plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='bead radius',limited=False,nsteps_list=np.arange(2,15,1),m_list=m_list,show=False,save=True,labelposition="outside",avoidLastStep=True)
    # Comparison between regular and limited, 10 steps
    plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='bead radius',bothLimitedAndNot=True,nsteps_list=[10],m_list=m_list,show=False,save=True,labelposition="outside",avoidLastStep=True,title = "Comparison between success rates for limited and regular walk,\n freely jointed chain with normal variated step length")

    # FREELY JOINTED NORMAL VARIATED, DON'T AVOID LAST BEAD
    m_list=np.arange(0.099,0.65,0.05)
    # Limited, comparsion of different numbers of steps
    plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='bead radius',limited=True,nsteps_list=np.arange(2,15,1),m_list=m_list,show=False,save=True,labelposition="outside",avoidLastStep=False)
    # Regular, comparsion of different numbers of steps
    plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='bead radius',limited=False,nsteps_list=np.arange(2,15,1),m_list=m_list,show=False,save=True,labelposition="outside",avoidLastStep=False)
    # Comparison between regular and limited, 10 steps
    plot_success_rate_vs_onXAxis([chainwalk_stepl_variations],onXAxis='bead radius',bothLimitedAndNot=True,nsteps_list=[10],m_list=m_list,show=False,save=True,labelposition="outside",avoidLastStep=False,title = "Comparison between success rates for limited and regular walk,\n freely jointed chain with normal variated step length,\n no collision between consecutive steps")

    ## COMPARISON BETWEEN FJ, FJ NORMVAR, LIMITED AND NOT,10 steps
    # Avoid last step
    plot_success_rate_vs_onXAxis([chainwalk,chainwalk_stepl_variations],onXAxis='bead radius',bothLimitedAndNot=True,nsteps_list=[10],m_list=np.arange(0,0.6,0.025),show=False,save=True,labelposition="outside",avoidLastStep=True)
    # Don't avoid last step
    plot_success_rate_vs_onXAxis([chainwalk,chainwalk_stepl_variations],onXAxis='bead radius',bothLimitedAndNot=True,nsteps_list=[10],m_list=np.arange(0,1.01,0.05),show=True,save=True,labelposition="outside",avoidLastStep=False)

    # plot_success_rate_vs_onXAxis([chainwalk],onXAxis='bead radius',limited=True,nsteps_list=[20],m_list=np.arange(0,2.01,0.025),show=True,save=True,labelposition="inside",avoidLastStep=False)
    #
    # # Check if success rate depends on r or just the relation between r and rho
    # plot_success_rate_vs_r([chainwalk],r_range=np.arange(1,30,1),M=[0.2,0.3,0.4],show=True,save=True,scale='linlin')

    ### END OF SUCCESS RATE PLOTS ###
    # ==========================================================================

    # #-------------------------# #
    # # # END-TO-END-DISTANCE # # #
    # #-------------------------# #

    # Find bead size that best corresponds to grid walk (not limited, self avoiding chain walk)
    # chainwalk.plot_multiple_bead_size_variations_holdon(nwalks=10000,nsteps=[5,7,9,11,13,15],types_of_walks=[1],data=[0],rhos=np.arange(0.099,0.5,0.05),plot_rms_grid=True,save=True,title="RMS of end-to-end distance for freely jointed chain\n and grid walker, 10000 walks",labels=[None],ylabel="RMS",boxes=False)

    # # Compare grid with FJ chain, not self avoiding
    # plot_multiple_end_to_end_distances([gridwalk,chainwalk],types_of_walks=[0],data=[0,1,2],nsteps_list=range(5,100,5),nwalks=5000,show=False,save=True)
    # # Compare grid with FJ chain, self avoiding
    # plot_multiple_end_to_end_distances([gridwalk,chainwalk],types_of_walks=[1],data=[0,1,2],nsteps_list=range(5,20,2),nwalks=5000,show=False,save=True)

    # # Compare self avoiding with limited self avoiding, grid
    # plot_multiple_end_to_end_distances([gridwalk],types_of_walks=[1,2],data=[0,1,2],nsteps_list=range(5,20,2),nwalks=5000,show=False,save=True)
    # # Compare limited self avoiding with forced self avoiding, grid
    # plot_multiple_end_to_end_distances([gridwalk],types_of_walks=[1,3],data=[0,1,2],nsteps_list=range(5,20,2),nwalks=5000,show=True,save=True)


    # Get rms for all grids of different number of steps (not important code)
    # plot_multiple_end_to_end_distances([gridwalk],types_of_walks=[1],data=[0],nsteps_list=range(1,20,1),nwalks=10000,show=False,save=True)

    # plot_multiple_end_to_end_distances([gridwalk,chainwalk,chainwalk_stepl_variations],types_of_walks=[0],data=[0],nwalks=1000,save=True)
    # plot_multiple_end_to_end_distances([gridwalk,chainwalk,chainwalk_stepl_variations],types_of_walks=[1],data=[1],nwalks=1000,save=True)
    # plot_multiple_end_to_end_distances([gridwalk,chainwalk,chainwalk_stepl_variations],types_of_walks=[2],data=[0],nwalks=1000,save=True)
    # plot_multiple_end_to_end_distances([gridwalk,chainwalk,chainwalk_stepl_variations],types_of_walks=[3],data=[0],nwalks=1000,save=True)
    # plot_multiple_end_to_end_distances([gridwalk,chainwalk,chainwalk_stepl_variations],types_of_walks=[0],data=[1],nwalks=1000,save=True)




    ### NORMPLOTS
    #GRID
    #gridwalk.normplot(nwalks=100,nsteps=100)
    #gridwalk.normplot(nwalks=100,nsteps=100,avoid=True,limited=True)

    #GRID STEPL VARIATIONS
    #grid_walker_stepl_variations.normplot(nwalks=100,nsteps=100)
    #grid_walker_stepl_variations.normplot(nwalks=100,nsteps=100,avoid=True,limited=True)

    #FJ
    #chainwalk.normplot(nwalks=100,nsteps=100)
    # chainwalk.normplot(nwalks=100, nsteps=15,avoid=True,limited=True)
    # chainwalk.normplot(nwalks=100, nsteps=15,avoid=True,limited=False)

    #FJ STEPL VARIATIONS
    #chainwalk_stepl_variations.normplot(nwalks=100,nsteps=100)
    #chainwalk_stepl_variations.normplot(nwalks=100, nsteps=15,avoid=True,limited=True)


    # gridwalk.plot_multiple_end_to_end_distances_holdon(show=False)
    # chainwalk.plot_multiple_end_to_end_distances_holdon()
main()
