import random as rnd
import numpy as np
import re
import matplotlib.pyplot as plt  # For 3D plot
from scipy.stats import norm  # histogram fitting
from scipy.stats import normaltest
from scipy import stats as scipy_stats
import math


class Walker():
    """Walked positions are stored in a list"""
    # Initiate list of visited points
    rho0 = 0.499  # size of first bead, the very least 1/2 step length.
    # Modify generate_rho() to manage method for generating the sizes of the other beads
    origin = [0, 0, 0, rho0]  # position and bead size are stored in list
    visited_points = [origin]
    length = 0
    my = 1
    sigma = 0.1

    # Define step length
    r = my

    # Store last walking direction
    last_direction = 0
    last_r = my

    def restart(self):
        """Resets list of visited points."""
        self.origin = self.origin[:3]
        self.origin.append(self.rho0)
        self.visited_points = [self.origin]
        self.last_direction = 0
        self.length = 0
        self.r = self.my


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
            theta_back = np.pi - self.last_direction[0]
            phi_back = np.pi + self.last_direction[1]

            # In local coordinate system with "going back" vector as z axis:
            # Generate new bead center
            alpha = np.arccos(
                (self.r ** 2 + self.last_r ** 2 - (2 * self.visited_points[-2][3]) ** 2) / (2 * self.r * self.last_r))
            z_rnd = rnd.uniform(-1, np.cos(alpha))
            phi = rnd.uniform(0, np.pi * 2)
            theta = np.arccos(z_rnd)

            # Translate bead center position into global coordinate system.
            # First: define the local coordinate system in terms of the global
            z_hat = [np.sin(theta_back) * np.cos(phi_back), np.sin(theta_back) * np.sin(phi_back),
                     np.cos(theta_back)]  # vector between the most recent two beads
            y_hat = [-np.sin(phi_back), np.cos(phi_back), 0]  # changed x->y
            # x_hat = [x_hat[i]* 1 / np.sin(theta_back) for i in range(3)]  # Orthogonal to z_hat
            x_hat = [-np.cos(theta_back) * np.cos(phi_back), -np.cos(theta_back) * np.sin(phi_back),
                     np.sin(theta_back)]  # Orthogonal to z_hat, x_hat #changed y->x
            #    y_hat = [y_hat[i]* 1 / np.sin(theta_back) for i in range(3)]  # Orthogonal to z_hat

            # print(sum([y_hat[i] * y_hat[i] for i in range(3)]))
            # Second: project the bead center position onto the global axes and translate it to origo
            current_pos_origo = [self.r * (
                        np.cos(theta) * z_hat[i] + np.sin(theta) * np.cos(phi) * y_hat[i] + np.sin(theta) * np.sin(
                    phi) * x_hat[i]) for i in range(3)]
            current_pos = [current_pos[i] + current_pos_origo[i] for i in range(3)]
            current_pos.append(rho)

            # vector_length = np.sqrt(sum([current_pos[i]**2 for i in range(3)]))
            # Define direction
            self.last_direction = [np.arctan(current_pos_origo[1] / current_pos_origo[0]),
                                   np.arccos(current_pos_origo[2] / self.r)]
            # self.last_direction = [np.arccos(current_pos[2]/ vector_length),np.arctan(current_pos[1] / current_pos[0])]
        else:
            z_rnd = rnd.uniform(-1, 1)
            phi = rnd.uniform(0, np.pi * 2)
            theta = np.arccos(z_rnd)

            self.last_direction = [theta, phi]
            # Update the coordinates
            current_pos[0] += self.r * np.sin(theta) * np.cos(phi)  # x
            current_pos[1] += self.r * np.sin(theta) * np.sin(phi)  # y
            current_pos[2] += z_rnd  # z
            current_pos[3] = rho
            ###---End of Suggestion---###
        # Update list of visited points
        self.visited_points.append(current_pos)
        self.length += self.r
        self.last_r = self.r
        return False

    def test_limited_option(self, limited=True, noptions=100):
        """Tests whether the limited-option actually limits options to > 2*rho from the neighbor of the most recent point"""

        def walk_options_help(limited=True):
            options = []
            current_pos = self.visited_points[-1][:]
            # Get bead size
            rho = self.rho0
            if limited == True and self.last_direction != 0:  # The last step

                # Define direction to walk back the same way
                theta_back = np.pi - self.last_direction[0]
                phi_back = np.pi + self.last_direction[1]
                # In local coordinate system with "going back" vector as z axis:
                # Generate new bead center
                self.alpha = np.arccos((self.r ** 2 + self.last_r ** 2 - (2 * self.visited_points[-2][3]) ** 2) / (
                            2 * self.r * self.last_r))
                # TESTING STEP: Generate many options for the position of the next bead. All options are stored in
                for i in range(noptions):
                    z_rnd = rnd.uniform(-1, np.cos(self.alpha))
                    phi = rnd.uniform(0, np.pi * 2)
                    theta = np.arccos(z_rnd)

                    # Translate bead center position into global coordinate system.
                    # First: define the local coordinate system in terms of the global
                    z_hat = [np.sin(theta_back) * np.cos(phi_back), np.sin(theta_back) * np.sin(phi_back),
                             np.cos(theta_back)]  # vector between the most recent two beads
                    y_hat = [-np.sin(phi_back), np.cos(phi_back), 0]  # Orthogonal
                    x_hat = [-np.cos(theta_back) * np.cos(phi_back), -np.cos(theta_back) * np.sin(phi_back),
                             np.sin(theta_back)]  # Orthogonal
                    # Second: project the bead center position onto the global axes and translate it to origo
                    option_direction = [self.r * (
                                np.cos(theta) * z_hat[i] + np.sin(theta) * np.cos(phi) * y_hat[i] + np.sin(
                            theta) * np.sin(phi) * x_hat[i]) for i in range(3)]
                    option = [current_pos[i] + option_direction[i] for i in range(3)]
                    options.append(option)
                # current_pos=option
                # current_pos.append(self.rho0)
            elif limited == False and self.last_direction != 0:  # The last step

                for i in range(noptions):
                    option = [0, 0, 0, self.rho0]
                    z_rnd = rnd.uniform(-1, 1)
                    phi = rnd.uniform(0, np.pi * 2)
                    theta = np.arccos(z_rnd)
                    option[0] = current_pos[0] + self.r * np.sin(theta) * np.cos(phi)  # x
                    option[1] = current_pos[1] + self.r * np.sin(theta) * np.sin(phi)  # y
                    option[2] = current_pos[2] + z_rnd  # z
                    options.append(option)
                # current_pos = option
                # current_pos.append(self.rho0)
            else:  # The first step
                z_rnd = rnd.uniform(-1, 1)
                phi = rnd.uniform(0, np.pi * 2)
                theta = np.arccos(z_rnd)

                self.last_direction = [theta, phi]
                # Update the coordinates
                current_pos[0] += self.r * np.sin(theta) * np.cos(phi)  # x
                current_pos[1] += self.r * np.sin(theta) * np.sin(phi)  # y
                current_pos[2] += z_rnd  # z
                current_pos[3] = rho
                # Update list of visited points
                self.visited_points.append(current_pos)
            self.length += self.r
            self.last_r = self.r
            return options

        # Generate the 3-pearl (2-steps) limited walk
        self.restart()
        for i in range(2):
            options = walk_options_help(limited)
        # Check to see the distance between the center of the first pearl and each option
        dist = []
        first_bead = self.visited_points[0]
        second_bead = self.visited_points[1]
        for opt in options:
            dist.append(np.sqrt(sum([(opt[i] - first_bead[i]) ** 2 for i in range(3)])))
        print("Distance between pearl 1 and options for pearl 3:\nmin: ", min(dist), ", 2*rho: ", 2 * self.rho0,
              ", max: ", max(dist), ", r:", self.r)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # I'm interested in seeing that the options span a full set of angles except for those that gives interception with the first pearl

        # Plot the vector originating from the option center, stretching rho0 length units in the direction toward the first bead pearl.
        # Result: For the non-limited option, we see that there are above vectors that penetrate the first pearl = interception, while those are not to be found in the limited option.
        # Result cont'd: At the same time, the limited option does not seem to restrict options unnecessarily much. As seen, there are green vectors that touch the first pearl's surface...
        # ...all over the side of the pearl. No angles around the first pearl is left out.

        # for opt in options:
        #     opt_firstBead = [first_bead[i] - opt[i] for i in range(3)]
        #     norm = np.sqrt(sum([opt_firstBead[i] ** 2 for i in range(3)]))
        #     opt_firstBead = [opt_firstBead[i] * (self.rho0 / norm) + opt[i] for i in range(3)]
        #     opt_firstBead_vectors = [[opt[i], opt_firstBead[i]] for i in
        #                              range(3)]  # vector from option to first bead with length rho0.
        #     #ax.plot(opt_firstBead_vectors[0], opt_firstBead_vectors[1], opt_firstBead_vectors[2], "green")
        #     #plt.suptitle("The green vectors originate from the center of the option,\n stretches a bead radius toward the center of the first pearl")
        #
        #     opt_secondBead_vectors = [[second_bead[i], opt[i]] for i in range(3)]  # vector from second bead to option (centerpoint)
        #     ax.plot(opt_secondBead_vectors[0], opt_secondBead_vectors[1], opt_secondBead_vectors[2], "orange")
        #     plt.suptitle("The orange vectors link the second pearl with the options")
        #     #plt.suptitle("The orange vectors link the second pearl with the options\nThe green vectors originate from the center of the option,\nstretches a bead radius toward the center of the first pearl")

        # Plot the options as spheres
        #         x = [i[0] for i in options]
        #         y = [i[1] for i in options]
        #         z = [i[2] for i in options]
        #         rho = self.rho0
        #         cmap = ["yellow", "green", "orange"]
        #         for i in range(len(x)):
        #             phi, theta = np.mgrid[0:2 * np.pi:200j, 0:np.pi:100j]
        #             x_sphere = x[i] + rho * np.cos(phi) * np.sin(theta)
        #             y_sphere = y[i] + rho * np.sin(phi) * np.sin(theta)
        #             z_sphere = z[i] + rho * np.cos(theta)
        #             ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color=cmap[2])
        #         plt.show()

        # ----- Plane intersection code -----#
        # Tried to project the option onto a plane between the first and second pearl (plane normal=step 1) but didn't work as I wanted it to. Copy-pasted
        # intersection function
        def isect_line_plane_v3(p0, p1, p_co, p_no, epsilon=10 ^ -6):
            """
            p0, p1: Define the line.
            p_co, p_no: define the plane:
                p_co Is a point on the plane (plane coordinate).
                p_no Is a normal vector defining the plane direction;
                     (does not need to be normalized).

            Return a Vector or None (when the intersection can't be found).
            """
            u = sub_v3v3(p1, p0)
            dot = dot_v3v3(p_no, u)
            if abs(dot) > epsilon:
                # The factor of the point between p0 -> p1 (0 - 1)
                # if 'fac' is between (0 - 1) the point intersects with the segment.
                # Otherwise:
                #  < 0.0: behind p0.
                #  > 1.0: infront of p1.
                w = sub_v3v3(p0, p_co)
                fac = -dot_v3v3(p_no, w) / dot
                u = mul_v3_fl(u, fac)
                return add_v3v3(p0, u)
            else:
                # The segment is parallel to plane.
                return None

        # ----------------------
        # generic math functions

        def add_v3v3(v0, v1):
            return [
                v0[0] + v1[0],
                v0[1] + v1[1],
                v0[2] + v1[2],
            ]

        def sub_v3v3(v0, v1):
            return (
                v0[0] - v1[0],
                v0[1] - v1[1],
                v0[2] - v1[2],
            )

        def dot_v3v3(v0, v1):
            return (
                    (v0[0] * v1[0]) +
                    (v0[1] * v1[1]) +
                    (v0[2] * v1[2])
            )

        def len_squared_v3(v0):
            return dot_v3v3(v0, v0)

        def mul_v3_fl(v0, f):
            return (
                v0[0] * f,
                v0[1] * f,
                v0[2] * f,
            )

        # ----- END Plane intersection code -----#
        plane_normal = [first_bead[i] - second_bead[i] for i in range(3)]
        norm = np.sqrt(sum([plane_normal[i] ** 2 for i in range(3)]))
        point_on_plane = [second_bead[i] + plane_normal[i] * self.r / (2 * norm) for i in range(3)]

        intersections = []
        for opt in options:
            opt_bead = [opt[i] - second_bead[i] for i in range(3)]
            if sum([opt_bead[i] * plane_normal[i]]) > 0:
                intersection = isect_line_plane_v3(second_bead, opt, point_on_plane, plane_normal)
                if intersection is not None and max([abs(intersection[i]) for i in range(3)]) < 2:
                    intersections.append(intersection)

        x = [i[0] for i in intersections]
        y = [i[1] for i in intersections]
        z = [i[2] for i in intersections]
        ax.plot(x, y, z, "ro")

        # Plot the first 2 pearls
        x = [i[0] for i in self.visited_points]
        y = [i[1] for i in self.visited_points]
        z = [i[2] for i in self.visited_points]
        rho = [i[3] for i in self.visited_points]

        ax.plot(x, y, z)
        cmap = ["yellow", "green", "orange"]
        for i in range(len(x)):  # if range(len(x)-1): just plot the first pearl
            phi, theta = np.mgrid[0:2 * np.pi:200j, 0:np.pi:100j]
            x_sphere = x[i] + rho[i] * np.cos(phi) * np.sin(theta)
            y_sphere = y[i] + rho[i] * np.sin(phi) * np.sin(theta)
            z_sphere = z[i] + rho[i] * np.cos(theta)
            ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color=cmap[i], alpha=0.2)

        plt.show()

    def test_avoid(self, avoidLastStep=True):
        """Test if latest site is already occupied - continuous case"""
        cp = self.visited_points[-1]  # Current position
        # print(cp)
        # The distance between successive sphere centres is self.r. Interception between any two spheres occurs if their centres are less apart than their diameter
        # if self.r < 2*self.rho:
        #     return True
        if avoidLastStep is True:
            last = -1
        elif avoidLastStep is False:
            last = -2

        for point in self.visited_points[:last]:
            r_centres = np.sqrt((point[0] - cp[0]) ** 2 + (point[1] - cp[1]) ** 2 + (point[2] - cp[2]) ** 2)
            if r_centres < point[3] + cp[3]:
                # Self-intercept - needs to restart the process
                return True
        return False
        # TODO Implement method for avoiding previous _paths_, not just previous positions.


def main():
    chainwalk = Freely_jointed_chain()

    chainwalk.test_limited_option(limited=True, noptions=10000)
    # chainwalk.test_limited_option(limited=False,noptions=10)


main()