import random as rnd
import numpy as np

class Grid_walker():
    """Walked positions are stored in three lists, one for each coordinate."""
    # Initial position
    x0 = 0
    y0 = 0
    z0 = 0

    # Define step length
    r = 1

    # Initiate the list storing the coordinates of the walk
    x, y, z = [x0], [y0], [z0]

    def walk_one_step():
        # Append the same coordinates as last step to coordinate list
        x.append(x[-1]), y.append(y[-1]), z.append(z[-1])
        # Update to new coordinates
        direction = rnd.randint(0,5)
        if direction == 0:
            x[-1] += 1
        elif direction == 1:
            x[-1] -= 1
        elif direction == 2:
            y[-1] += 1
        elif direction == 3:
            y[-1] -= 1
        elif direction == 4:
            z[-1] += 1
        elif direction == 5:
            z[-1] -= 1

class Freely_jointed_chain():
    """Walked positions are stored in three lists, one for each coordinate."""
    # Initial position
    x0 = 0
    y0 = 0
    z0 = 0

    # Define step length
    r = 1

    # Initiate the list storing the coordinates of the walk
    x, y, z = [x0], [y0], [z0]

    def walk_one_step():
        # Append the same coordinates as last step to coordinate list
        x.append(x[-1]), y.append(y[-1]), z.append(z[-1])
        # Get walking direction
        theta = rnd.uniform(0,np.pi)
        phi = rnd.uniform(0,2*np.pi)
        # Update the coordinates
        x[-1] += r*np.sin(theta)*np.cos(phi)
        y[-1] += r*np.sin(theta)*np.sin(phi)
        z[-1] += r*np.cos(theta)
