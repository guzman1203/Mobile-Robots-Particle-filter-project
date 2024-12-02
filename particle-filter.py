import math
import time
import numpy as np
from scipy.stats import norm
from fairis_tools.my_robot import MyRobot







# Initialize robot and sensors
robot = MyRobot()

GRID_SIZE = 5  # 5x5 grid
CELL_SIZE = 1.0  # Each grid cell is 1x1 meter
MAX_CELLS = 20

# Max and Min motor velocities
MAX_ROTATIONAL = 4.0
MAX_VELOCITY = robot.max_motor_velocity

# Cardinal Orientations
EAST = 0
NORTH = 90
WEST = 180
SOUTH = 270
CARDINALS = [EAST, NORTH, WEST, SOUTH]

# tolerance on PIDs
THRESHOLD_FROM_WALL = 1
ORIENTATION_TOLERANCE = math.radians(0.1)
MOVEMENT_TOLERANCE = 0.001

# robot camera properties
CAMERA_WIDTH = robot.rgb_camera.getWidth()
CAMERA_CENTER_X = CAMERA_WIDTH / 2

# Load maze environment
robot.load_environment('../../worlds/Fall24/maze8.xml')
# Get the time step 
timestep = int(robot.getBasicTimeStep())
# Move robot to a random staring position 
robot.move_to_start()
 
# Known wall configuration
WALL_CONFIG = [
    ['O', 'W', 'W', 'O'], ['O', 'W', 'O', 'W'], ['O', 'W', 'O', 'W'], ['O', 'W', 'O', 'W'], ['W', 'W', 'O', 'O'],
    ['W', 'O', 'W', 'O'], ['O', 'W', 'W', 'W'], ['O', 'W', 'O', 'W'], ['O', 'W', 'O', 'W'], ['W', 'O', 'O', 'W'],
    ['O', 'O', 'W', 'O'], ['O', 'W', 'O', 'W'], ['O', 'W', 'O', 'W'], ['O', 'W', 'O', 'W'], ['W', 'W', 'O', 'O'],
    ['O', 'O', 'W', 'O'], ['O', 'W', 'O', 'W'], ['O', 'W', 'O', 'W'], ['O', 'W', 'O', 'W'], ['W', 'O', 'O', 'O'],
    ['O', 'O', 'W', 'W'], ['O', 'W', 'O', 'W'], ['O', 'W', 'O', 'W'], ['O', 'W', 'O', 'W'], ['W', 'O', 'O', 'W'],
]

# Navigation through the maze
    # the robot needs to traverse each cell one by one to allow for a sensor readings
    # new cells can be detected VIA
        # moving 1 meter
    # the robot needs to tell the difference between 15 cells
    # the robot would be set up with a dead end before the 15 different cells
        # uniques = cells with unique patterns
        # available tools
            # compass
            # encoders
            # lidar
        # 
        
        # ideas
            # finding a unique orientates location
            # keeping track of potential paths using an array of potential cells, where each cell represents the head of the path
                # we keep track of the movements the robot makes
                # for every new cell traveled to, 
                    # we apply the movement to our path cells
                    # we measure the cells walls
                    # we cull impossible paths

        # what we need
            # to keep track of the movements we take
            # to tell when we reached a dead end
            # to tell when we backtrack on previously tracked ground
            
            # keep a group of path head cells
            # be able to apply transitions to path head cells, leading to path simulation
            # cull mismatching path heads given current sensor model readings
            
            # 
            
        # what we have
            # we have the ability to get a list of cells that we have the most probability of residing in
            # get a list of cardinal sensor model readings
            # to move around the maze
         
        
        # movement
            # going around corners logic
                # 
            
        # movement tracking
            # after every PID move forward, append an F to the string
            # after every setOrientation(WNES), append the corresponding cardinal to the string
            
            # could set the string to repeating letter format
            
            # W-F-F-F-F-N-F-E-F-F-F-F-N-F-W-F-F-F-F
            
            
        # starting a graph
            # use a graph and a visited coordinate set to keep track of visited cells
            # when the robot starts, it also starts a graph with 1 node with a coordinate of 0,0, also the visited set has 0,0
            # as the robot enters a new cell, it adds to the graph and the visited set based on the direction it took
            
            # graph is implemented using a dictionary with coordinates as keys for easy lookup
            # the graph initializes at a cell node with 0,0 and branches out with new cell nodes from there
            
            # structure
                # graph
                    # a global changable dictionary to store coordinate keys nodes value pairs for easy lookup
                    # a get_node(coordinate) function that returns a node given a coordinate
                
                # node
                    # key: coorinates
                    # neighbor list whitch point to other nodes
                        # the edges can be infered          
            
            # initialization
                # the first node is at 0,0 for the graph 
                # we add leafnodes/neighbors corresponding to open walls in the cardinal directions
                # we can calc the neighbors coordinates via:
                    # east  = y   , x+1
                    # north = y+1 , x
                    # west  = y   , x-1
                    # south = y-1 , x
                # keep a head node with the next positionition of the robot in the graph
                     
            # for every new cell
                # 
                # get the leafnode from the graph
                # 
                # the graph calculates 
                
                
        # visited set implementation pseudocode
            # summary
            # use a visited set of coordinate keys starting at 0,0 to keep track of where the robot has been
            # assume movement logic will traverse all cells given this visited set 
            # will use keys in r-c format
            
            # initialization
            # start the visited set with the element 0-0
            # keep a headcell with the current robots coordinates
            
            # for every new cell
                # use the cardinal taken to calc the new_cell's position coordinates
                    # we know the cardinal from after moving a cell and calling get_current_facing_cardinal()
                    # moved East  :  y   , x+1  
                    # moved North :  y+1 , x
                    # moved West  :  y   , x-1
                    # moved south :  y-1 , x-1
                # add a new key to the visited set with those coordinates
        
        # movement logic
            # summary - based on start location
            # the idea is to start with keeping a start coordinate, following the left wall, and a need_new_start flag and left_wall_following flag 
            # at every cell determine what is the coordinate you are traversing to next
            # ask if it is that cell is the start coordinate
                # if it is
                    # switch wall followings
                    # and toggle need_new_start
                # if it has NOT
                    # keep following the same wall
                    # if need_new_start
                        # toggle need_new_start
                        # start = these coordinates

            # summary - based on visited locations
            # the idea is to start with following the left wall and ignore_visited flag 
            # at every cell determine what is the coordinate you are traversing to next
            # ask if it is that cell has been visited
                # if it is and ignore_visited == False
                    # switch wall followings
                    # ignore_visited = True
                # if it has NOT
                    # ignore_visited = False
                    # keep following the same wall
                    

# COMPLETE
    # movement
    # wall detecting in correct orientation
    # cell_probabilities and get_highest_prob
    
# NOT COMPLETE
    # perfecting movement
    # traversal logic
    # tracking traversal
        # traversal string
        # or graph
        # or visited coordinates



# Function to orient the robot to a target orientation
def set_orientation_pid(target_orientation):
    """
    Orient the robot to a specified orientation (in radians) using a PID controller,
    with motor velocity limits of ±4, and using robot.step(timestep) for Webots compatibility.

    Args:
        target_orientation_rad (float): The desired orientation in radians (0 to 2π).
        timestep (int): The simulation timestep for robot.step.

    Returns:
        bool: True if the robot successfully orients itself within a tolerance, False otherwise.
    """
    
    # Convert target orientation from degrees to radians
    target_orientation_rad = math.radians(target_orientation)
    
    # PID control parameters
    Kp = 8.0   # Proportional gain
    Ki = 0.01  # Integral gain
    Kd = 2.0   # Derivative gain

    # Error values for PID
    prev_error = 0.0
    integral = 0.0
    dt = 0.032


    while robot.step(timestep) != -1:
            
        # Get the robot's current orientation
        current_orientation = math.radians(robot.get_compass_reading())  # In radians
        #print(f"current_orientation={current_orientation:.4f}m")

        # Calculate the error in orientation
        error = target_orientation_rad - current_orientation
        #print(f"error={error:.4f}m")

        # If the error is within tolerance, stop the robot and return success
        if abs(error) <= ORIENTATION_TOLERANCE:
            print(f"    Robot oriented to {math.degrees(target_orientation_rad):.4f} degrees.")
            robot.set_left_motors_velocity(0)
            robot.set_right_motors_velocity(0)
            return True

        # PID calculations
        proportional = Kp * error
        integral += error * Ki * dt
        derivative = Kd * (error - prev_error) / dt

        # Compute control signal
        control_signal = proportional + integral + derivative
        prev_error = error
        
        # Reduce motor speed dynamically as error decreases
        if abs(error) < math.radians(10):  
            control_signal *= 0.5

        # Limit the control signal to the maximum motor velocity range [-4, 4]
        control_signal = max(min(control_signal, MAX_ROTATIONAL), -MAX_ROTATIONAL)
        #print(f"control_signal={control_signal:.4f}m")

        # Set motor velocities based on control signal
        robot.set_left_motors_velocity(-control_signal)
        robot.set_right_motors_velocity(control_signal)
        
        # print(f"     current orienation {math.degrees(current_orientation):.4f}     control signal: {control_signal:.2f}\n")

def set_orientation_east():
    print(f"Orienting to the EAST:{EAST}")
    set_orientation_pid(EAST)

def set_orientation_north():
    print(f"Orienting to the NORTH:{NORTH}")
    set_orientation_pid(NORTH)

def set_orientation_west():
    print(f"Orienting to the WEST:{WEST}")
    set_orientation_pid(WEST)

def set_orientation_south():
    print(f"Orienting to the SOUTH:{SOUTH}")
    set_orientation_pid(SOUTH)

def normalize_degree(angle):
    """Normalize an angle to be within the range [-π, π]."""
    return (math.radians(angle) + math.pi) % (2 * math.pi) - math.pi

def set_orientation_cardinal_reverse():
    print(f"Reversing orientation")
    current = get_current_facing_cardinal()
    opposites = {EAST: WEST, NORTH: SOUTH, 
                 WEST: EAST, SOUTH: NORTH}
    set_orientation_pid(opposites[current])

def set_orientation_turn_left_90():
    left_turns = {
        EAST  : NORTH,
        NORTH : WEST,
        WEST  : SOUTH,
        SOUTH : EAST
    }
    facing_cardinal = get_current_facing_cardinal()
    set_orientation_pid(left_turns[facing_cardinal])
    
def set_orientation_turn_right_90():
    right_turns = {
        EAST  : SOUTH,
        NORTH : EAST,
        WEST  : NORTH,
        SOUTH : WEST
    }
    facing_cardinal = get_current_facing_cardinal()
    set_orientation_pid(right_turns[facing_cardinal])
    
def get_current_facing_cardinal():
    current = normalize_degree(robot.get_compass_reading())
    #print(f"     normalized current: {math.degrees(current)}")
    #print(f"     list of cardinals: E:{math.degrees(normalize_degree(EAST))}   N:{math.degrees(normalize_degree(NORTH))}   W:{math.degrees(normalize_degree(WEST))}   -W:{math.degrees(normalize_degree(-WEST))}   S:{math.degrees(normalize_degree(SOUTH))}")    
    
    cardinals = {EAST : normalize_degree(EAST), 
                 NORTH: normalize_degree(NORTH), 
                 -WEST: normalize_degree(WEST),
                 WEST : -normalize_degree(WEST),
                 SOUTH: normalize_degree(SOUTH)}
                 
    closest_cardinal = min(cardinals, key=lambda cardinal: abs(cardinals[cardinal] - current))
    #print(f"     closest cardinal: {closest_cardinal}")
    
    return abs(closest_cardinal)

# Fuction to move the robot a target distance forward
def move_forward_pid(target_distance, start_timestep):
    """
    Moves the robot forward by a specified target distance (in meters) using a PID controller
    based on the front-left motor encoder readings.

    Args:
        target_distance (float): The target distance to move forward, in meters.
        timestep (int): The simulation timestep for robot.step.
    """
    # PID control parameters
    Kp = 15.0  # Proportional gain
    Ki = 0.05  # Integral gain
    Kd = 0.5   # Derivative gain

    # Encoder readings
    initial_encoder = robot.get_front_left_motor_encoder_reading()

    # Error values for PID
    prev_error = 0
    integral = 0
    dt = 0.032

    # print(f"     Starting PID-based forward movement: Target Distance={target_distance:.2f} meters")

    while robot.step(robot.timestep) != -1:
            
        # Get the current encoder reading
        current_encoder = robot.get_front_left_motor_encoder_reading()

        # Calculate the distance traveled based on the encoder
        radians_traveled = current_encoder - initial_encoder
        distance_traveled = robot.wheel_radius * radians_traveled

        # Calculate the error
        error = target_distance - distance_traveled

        # If the robot has traveled the target distance within the tolerance, stop
        if abs(error) <= MOVEMENT_TOLERANCE:
            robot.set_left_motors_velocity(0)
            robot.set_right_motors_velocity(0)
            # print(f"     Movement completed: Traveled Distance={distance_traveled:.2f} meters")
            return True

        # PID calculations
        proportional = Kp * error
        integral += error * dt  # Convert timestep to seconds
        derivative = Kd * (error - prev_error) / dt
        prev_error = error

        # Compute control signal
        control_signal = proportional + integral + derivative

        # Limit control signal to maximum motor velocity
        control_signal = max(min(control_signal, MAX_VELOCITY), -MAX_VELOCITY)

        # Set motor velocities to move forward
        robot.set_left_motors_velocity(control_signal)
        robot.set_right_motors_velocity(control_signal)


        # Debugging output
        # print(f"     Encoder Reading: {current_encoder:.2f}, Distance Traveled: {distance_traveled:.2f}, "
        #      f"     Error: {error:.2f}, Control Signal: {control_signal:.2f}")
    #print(f"     Current Servomotors Velocity: {control_signal :.2f} m/s")
    
def move_forward_one_cell():
    print(f"Moving forward 1 cell")
    move_forward_pid(CELL_SIZE, robot.timestep) 

def move_forward_x_cells(x):
    for i in range(x):
        move_forward_one_cell()

def get_front_lidar_reading():
    return robot.get_lidar_range_image()[400]
    
def get_left_lidar_reading():
    return robot.get_lidar_range_image()[200]

def get_right_lidar_reading():
    return robot.get_lidar_range_image()[600]

def get_rear_lidar_reading():
    return robot.get_lidar_range_image()[0]

def get_oriented_lidar_readings():
    current_cardinal = get_current_facing_cardinal()
    result = False
    
    lidar_readings = {
        EAST: {
            EAST : get_front_lidar_reading(),
            NORTH: get_left_lidar_reading(),
            WEST : get_rear_lidar_reading(),
            SOUTH: get_right_lidar_reading(),
        },
        NORTH: {
            EAST : get_right_lidar_reading(),
            NORTH: get_front_lidar_reading(),
            WEST : get_left_lidar_reading(),
            SOUTH: get_rear_lidar_reading(),
        },
        WEST: {
            EAST : get_rear_lidar_reading(),
            NORTH: get_right_lidar_reading(),
            WEST : get_front_lidar_reading(),
            SOUTH: get_left_lidar_reading(),
        },
        SOUTH: {
            EAST : get_left_lidar_reading(),
            NORTH: get_rear_lidar_reading(),
            WEST : get_right_lidar_reading(),
            SOUTH: get_front_lidar_reading(),
        }
    }
    
    # Ensure that the current cardinal is valid
    if current_cardinal not in lidar_readings:
        raise ValueError(f"Invalid direction: {current_cardinal}.")
        
    return lidar_readings[current_cardinal]

def is_dead_end():
    return get_front_lidar_reading() < THRESHOLD_FROM_WALL and \
           get_left_lidar_reading() < THRESHOLD_FROM_WALL and \
           get_right_lidar_reading() < THRESHOLD_FROM_WALL
    
def get_walls(direction):
    lidar_readings = get_oriented_lidar_readings()
    walls_dict = {cardinal: lidar_readings[cardinal] < THRESHOLD_FROM_WALL
                  for cardinal in lidar_readings}
    return walls_dict

def is_front_wall():
    return get_front_lidar_reading() < THRESHOLD_FROM_WALL

def is_left_wall():
    return get_left_lidar_reading() < THRESHOLD_FROM_WALL
    
def is_right_wall():
    return get_right_lidar_reading() < THRESHOLD_FROM_WALL
    
def is_rear_wall():
    return get_rear_lidar_reading() < THRESHOLD_FROM_WALL
    
def print_walls(direction):
    print(f"Direction of reading: {direction}")
    walls = get_walls(direction)
    
    # Top and bottom borders, adjust depending on walls
    top_wall = " -------"  # Default for no wall
    if walls[NORTH]:  # North wall detected
        top_wall = " XXXXXXX"
    
    bottom_wall = " -------"  # Default for no wall
    if walls[SOUTH]:  # South wall detected
        bottom_wall = " XXXXXXX"
    
    # Side walls (left and right)
    west_wall = 'X' if walls[WEST] else '|'
    east_wall = 'X' if walls[EAST] else '|'
    
    # Now, build the string for the box:
    print(top_wall)  # Top border
    
    # Middle parts of the box (left and right walls)
    for _ in range(3):  # Print three middle lines with possible walls
        print(f"{west_wall}       {east_wall}")
    
    print(bottom_wall)  # Bottom border

# Sensor model
def calculate_cell_probabilities(wall_readings):
    probabilities = []
    for cell_idx, cell_config in enumerate(WALL_CONFIG):
        #print(f"     Cell #: {cell_idx+1}    config:{cell_config}")
        cell_prob = 1.0
        for i, cardinal in enumerate(CARDINALS) :
            s = 1 if cell_config[i] == 'W' else 0
            z = 1 if wall_readings[cardinal] is True else 0   
            #print(f"     cardinal: {cardinal}")
            #print(f"     z = {z} | s = {1}")
            if z == 1 and s == 1:
                cell_prob *= 0.8
            elif z == 1 and s == 0:
                cell_prob *= 0.4
            elif z == 0 and s == 1:
                cell_prob *= 0.2
            elif z == 0 and s == 0:
                cell_prob *= 0.6
        probabilities.append(round(cell_prob,5))
    #print(probabilities)
    
    # normalize probabilities
    sum_of_probs = sum(probabilities)
    norm_probs = [prob/sum_of_probs for prob in probabilities]
    
    return norm_probs

# Gather highest weighted cell probabilities
def get_highest_probability_cells(probabilities):
    max_prob = max(probabilities)
    highest_prob_cells = [i + 1 for i, prob in enumerate(probabilities) if prob == max_prob]
    return highest_prob_cells

def print_probabilities(probabilities):

    '''TODO: Prin in a grid'''

    for cell, prob in enumerate(probabilities, start=1):
        print(f"Cell {cell}: Probability = {prob:.3f}")

def get_next_position(curr_position):
    """
    Calculate the next position to move to based on the next positionition.

    Args:
        next_positionition (tuple): The robot's current grid position as (r, c).

    Returns:
        tuple: The next grid position as (r, c).
    """
    r, c = curr_position
    cardinal_direction = get_current_facing_cardinal()
    result = None

    if cardinal_direction == EAST    : result = (r, c + 1)
    elif cardinal_direction == NORTH : result = (r + 1, c)
    elif cardinal_direction == WEST  : result = (r, c - 1)
    elif cardinal_direction == SOUTH : result = (r - 1, c)

    return result

def get_key_from_position(position):
    return f"{position[0]}-{position[1]}"

def tracking_path_logic(visited, ignore_visited, nxt_position):
    # summary - based on visited locations
    # the idea is to start with following the left wall and ignore_visited flag 
    # at every cell determine what is the coordinate you are traversing to next
    # ask if it is that cell has been visited
        # if it is and ignore_visited == False
            # switch wall followings
            # ignore_visited = True
        # if it has NOT
            # ignore_visited = False
            # keep following the same wall
    
    wall_following_changes = False
    
    if get_key_from_position(nxt_position) in visited and not ignore_visited:
        wall_following_changes = True
        ignore_visited = True
    else:
        ignore_visited = False

    return wall_following_changes, ignore_visited

# Robot navigates through the maze function
def navigate_maze(robot):
    
    path_head = (0,0)
    visited = set([get_key_from_position(path_head)])
    ignore_visited = False
    following_left_wall = True

    while robot.step(robot.timestep) != -1 and len(visited) < MAX_CELLS:
        
        wall_following_changes = False
        
        print(f"Current path head is {path_head}")
        print(f"Following the {'left' if following_left_wall else 'right'} wall")
        
        # if dead end reverse
        if is_dead_end():
            print("Deadend detected. Reversing.")
            set_orientation_cardinal_reverse()
        
        elif following_left_wall:
            if not is_left_wall():
                print("Left wall not detected. Turning left.")
                set_orientation_turn_left_90()
            elif is_front_wall():
                print("Front wall detected. Turning right.")
                set_orientation_turn_right_90()
        
        elif not following_left_wall:
            if not is_right_wall():
                print("Right wall not detected. Turning right.")
                set_orientation_turn_right_90()
            elif is_front_wall():
                print("Front wall detected. Turning left.")
                set_orientation_turn_left_90()
        
        
        # before the robot moves forward
        # read out probable cells
        walls = get_walls(get_current_facing_cardinal())
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        probable_cells = get_highest_probability_cells(probs)
        print(f"Predicted cells are: {probable_cells}")
        
        # calc next position
        next_position = get_next_position(path_head)
        
        # determine path logic
        wall_following_changes, ignore_visited = tracking_path_logic(visited, ignore_visited, next_position)
        
        # toggle which wall the robot is following if needed
        if wall_following_changes:
            print(f"Visited cell ahead. Changing wall following from {'left' if following_left_wall else 'right'} to {'left' if not following_left_wall else 'right'}")
            following_left_wall = not wall_following_changes
        
        # move forward
        move_forward_one_cell()
        
        # if not visited, add next position to visited
        coord_key = get_key_from_position(next_position)
        if coord_key not in visited:
            print("New cell detected. Adding to visited cells.")
            visited.add(coord_key)
        
        # update path head
        path_head = next_position
        

    print("Robot visited 15 different cells!")
    
def test_components(robot):
    
    while robot.step(robot.timestep) != -1:
        
        ''' # testing movement
        set_orientation_west()
        move_forward_x_cells(4)
        
        set_orientation_north()
        move_forward_x_cells(1)
        
        set_orientation_east()
        move_forward_x_cells(4)
        
        set_orientation_north()
        move_forward_x_cells(1)
        
        set_orientation_west()
        move_forward_x_cells(4)
        
        set_orientation_north()
        move_forward_x_cells(2)
        
        set_orientation_east()
        move_forward_x_cells(4)
        
        set_orientation_south()
        move_forward_x_cells(1)
        
        set_orientation_west()
        move_forward_x_cells(3)
        '''
        ''' # testing reverse orientation
        set_orientation_east()
        print(f"    real:{robot.get_compass_reading()} cardinal:{get_current_facing_cardinal()}")
        set_orientation_cardinal_reverse()
        print(f"    real:{robot.get_compass_reading()} cardinal:{get_current_facing_cardinal()}")
    
        set_orientation_west()
        print(f"    real:{robot.get_compass_reading()} cardinal:{get_current_facing_cardinal()}")
        set_orientation_cardinal_reverse()
        print(f"    real:{robot.get_compass_reading()} cardinal:{get_current_facing_cardinal()}")
    
    
        set_orientation_north()
        print(f"    real:{robot.get_compass_reading()} cardinal:{get_current_facing_cardinal()}")
        set_orientation_cardinal_reverse()
        print(f"    real:{robot.get_compass_reading()} cardinal:{get_current_facing_cardinal()}")
   
        
        set_orientation_south()
        print(f"    real:{robot.get_compass_reading()} cardinal:{get_current_facing_cardinal()}")
        set_orientation_cardinal_reverse()
        print(f"    real:{robot.get_compass_reading()} cardinal:{get_current_facing_cardinal()}")
        '''
        '''# test wall detection
        set_orientation_west()
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        
        set_orientation_north()
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        
        set_orientation_east()
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())

        set_orientation_north()
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        
        set_orientation_west()
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
 
        set_orientation_north()
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        
        set_orientation_east()
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        
        set_orientation_south()
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        
        set_orientation_west()
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        move_forward_x_cells(1)
        print_walls(get_current_facing_cardinal())
        '''
        '''# test calculate_cell_probabilities and get_highest_probability_cells
        set_orientation_west()
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        
        set_orientation_north()
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        
        set_orientation_east()
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        
        set_orientation_north()
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        
        set_orientation_west()
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        
        set_orientation_north()
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        
        set_orientation_east()
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        
        set_orientation_south()
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        
        set_orientation_west()
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        move_forward_x_cells(1)
        probs = calculate_cell_probabilities(get_walls(get_current_facing_cardinal()))
        print(get_highest_probability_cells(probs))
        '''
        '''# test getting next position
        set_orientation_west()
        next_position = get_next_position((0,0))
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        

        set_orientation_north()
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        

        set_orientation_east()
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        
        
        set_orientation_north()
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        
        
        set_orientation_west()
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        
        set_orientation_north()
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        
        
        set_orientation_east()
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        
        
        set_orientation_south()
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        
        
        set_orientation_west()
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        next_position = get_next_position(next_position)
        print(f"next position: {next_position}")
        move_forward_x_cells(1)
        '''

        
        # test path tracking and navigation
        
    '''
    visited_cells = set()
    for _ in range(15):  # Visit at least 15 cells
        readings = get_sensor_readings(robot)  # sensor readings in cardinal directions
        probabilities = calculate_cell_probabilities(readings)
        
        print_probabilities(probabilities)
        
        highest_prob_cells = get_highest_probability_cells(probabilities)
        print(f"Predicted cells with the highest probability: {highest_prob_cells}")
        
        
        current_cell = get_current_cell(robot)  # Simulated function to determine current cell ########
        visited_cells.add(current_cell)
        
        
        move_to_next_cell(robot, visited_cells)  # Move to a new cell
    '''
    

# Main program
def main():
    # Main loop: Perform simulation Iterations until Webots stops the controller
    start_time = time.time()

    navigate_maze(robot)
    
    # Calculate total travel time
    end_time = time.time()
    travel_time = end_time - start_time
    print(f"Total travel time: {travel_time:.2f} seconds\n ")

if __name__ == "__main__":
    main()