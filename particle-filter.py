import math
import time
from fairis_tools.my_robot import MyRobot

# Initialize robot and sensors
robot = MyRobot()

# Robot constants
MAX_ROTATIONAL = 4.0
MAX_VELOCITY = robot.max_motor_velocity

# Camera constants
CAMERA_WIDTH = robot.rgb_camera.getWidth()
CAMERA_CENTER_X = CAMERA_WIDTH / 2

# Environment
GRID_SIZE = 5  # 5x5 grid
CELL_SIZE = 1.0  # cell size is 1 meter by 1 meter

# Cardinal orientations
EAST = 0
NORTH = 90
WEST = 180
SOUTH = 270
CARDINALS = [EAST, NORTH, WEST, SOUTH]

# PID controller tolerances
THRESHOLD_FROM_WALL = 1
ORIENTATION_TOLERANCE = math.radians(0.1)
MOVEMENT_TOLERANCE = 0.001

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

# Function to orient the robot to a target orientation
def pid_set_orientation(target):
    """
    Orient the robot to a specified orientation using a PID controller,
    using the robot's compass utility and accomplished via the robots motor velocities.

    Args:
        target_orientation_rad (float): The desired orientation in degrees.
    """
    
    # Check to change offet to compass readings for orientation to take the shortest path
    compass_offset = 0
    if abs(target - robot.get_compass_reading()) > 180:
        target = (target + 180) % 360
        compass_offset = 180
        
    # Convert target orientation from degrees to radians
    target_rads = math.radians(target)
    
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
        current_orientation = math.radians((robot.get_compass_reading() + compass_offset) % 360)  # In radians
        #print(f"current_orientation={current_orientation:.4f}m")

        # Calculate the error in orientation
        error = target_rads - current_orientation
        #print(f"error={error:.4f}m")

        # If the error is within tolerance, stop the robot and return success
        if abs(error) <= ORIENTATION_TOLERANCE:
            print(f"    Robot oriented to {math.degrees(target_rads):.4f} degrees.")
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
    pid_set_orientation(EAST)

def set_orientation_north():
    print(f"Orienting to the NORTH:{NORTH}")
    pid_set_orientation(NORTH)

def set_orientation_west():
    print(f"Orienting to the WEST:{WEST}")
    pid_set_orientation(WEST)

def set_orientation_south():
    print(f"Orienting to the SOUTH:{SOUTH}")
    pid_set_orientation(SOUTH)

def normalize_degree(angle):
    """Normalize an angle to be within the range [-π, π]."""
    return (math.radians(angle) + math.pi) % (2 * math.pi) - math.pi

def set_orientation_cardinal_reverse():
    print(f"Reversing orientation")
    current = get_current_facing_cardinal()
    opposites = {EAST: WEST, NORTH: SOUTH, 
                 WEST: EAST, SOUTH: NORTH}
    pid_set_orientation(opposites[current])

def set_orientation_turn_left_90():
    left_turns = {
        EAST  : NORTH,
        NORTH : WEST,
        WEST  : SOUTH,
        SOUTH : EAST
    }
    facing_cardinal = get_current_facing_cardinal()
    pid_set_orientation(left_turns[facing_cardinal])
    
def set_orientation_turn_right_90():
    right_turns = {
        EAST  : SOUTH,
        NORTH : EAST,
        WEST  : NORTH,
        SOUTH : WEST
    }
    facing_cardinal = get_current_facing_cardinal()
    pid_set_orientation(right_turns[facing_cardinal])
    
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
def pid_move_forward(target_distance, start_timestep):
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
    pid_move_forward(CELL_SIZE, robot.timestep) 

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

def get_walls(direction):
    lidar_readings = get_oriented_lidar_readings()
    walls_dict = {cardinal: lidar_readings[cardinal] < THRESHOLD_FROM_WALL
                  for cardinal in lidar_readings}
    return walls_dict

def is_dead_end():
    return get_front_lidar_reading() < THRESHOLD_FROM_WALL and \
           get_left_lidar_reading() < THRESHOLD_FROM_WALL and \
           get_right_lidar_reading() < THRESHOLD_FROM_WALL
    
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

    '''TODO: Print in a grid'''

    for cell, prob in enumerate(probabilities, start=1):
        print(f"Cell {cell}: Probability = {prob:.3f}")

def get_next_position(current_position):
    """
    Calculate the next position to move to based on the next positionition.

    Args:
        next_positionition (tuple): The robot's current grid position as (r, c).

    Returns:
        tuple: The next grid position as (r, c).
    """
    r, c = current_position
    cardinal_direction = get_current_facing_cardinal()
    result = None

    if cardinal_direction == EAST    : result = (r, c + 1)
    elif cardinal_direction == NORTH : result = (r + 1, c)
    elif cardinal_direction == WEST  : result = (r, c - 1)
    elif cardinal_direction == SOUTH : result = (r - 1, c)

    return result

def get_key_from_position(position):
    return f"{position[0]}-{position[1]}"

def path_logic(visited, ignore_visited, nxt_position):
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

    while robot.step(robot.timestep) != -1 and len(visited) < 15:
        
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
        wall_following_changes, ignore_visited = path_logic(visited, ignore_visited, next_position)
        
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
        

    print("Robot visited {MAX_CELLS} different cells!")
    
# Main program
def main():
    # Main loop: Perform simulation Iterations until Webots stops the controller
    start_time = time.time()

    navigate_maze(robot)
    # while robot.timestep != -1:
    #     move_forward_one_cell()
    #     set_orientation_north()
    #     set_orientation_east()
    #     set_orientation_west()
    #     set_orientation_south()


    # Calculate total travel time
    end_time = time.time()
    travel_time = end_time - start_time
    print(f"Total travel time: {travel_time:.2f} seconds\n ")

if __name__ == "__main__":
    main()