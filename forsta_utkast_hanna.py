import random as rnd

class Grid_walker():
    # Initial position
    x0 = 0
    y0 = 0
    z0 = 0

    # Initiate the list storing the coordinates of the walk
    x, y, z = [x0], [y0], [z0]

    def walk_one_step():
        x.append(x[-1]), y.append(y[-1]), z.append(z[-1])   # Append the same coordinates as last step to coordinate list
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
    initial_position = (0,0)

    def walk_one_step():
        pass
