import random as rnd
import numpy as np
limited=True
saw = [[0,0,0],[1,0,0],[1,1,0]]

# Choose an end point at random
choice=rnd.choice([0,-1])
print('one end',choice)
print('other',choice*(-1)-1)
# Remove this end
print('original',saw)
saw.pop(choice)  # One end
print('reduced',saw)
# Add a step on the other end identified as "current_pos"
current_pos = list(saw[choice * (-1) - 1])  # The other end
prev_pos = saw[choice * (-1) + (choice + 1) * (-2)]  # om choice i -1: prev_pos index= 1; om choice i 0: prev_pos index= -2
# Get walking direction
possible_directions = [-3, -2, -1, 1, 2, 3]
direction = rnd.choice(possible_directions)
# Get the last direction
for i in range(3):
    step = current_pos[i] - prev_pos[i]  # Step is either +1 or -1
    if abs(step) > 0:
        last_direction = (i + 1) * step
        break
if limited == True:  # Should always be True
    while direction == -last_direction:
        direction = rnd.choice(possible_directions)
print(last_direction,direction)
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
if choice == -1:
    saw.insert(0, current_pos)
else:
    saw.append(current_pos)
print(saw)