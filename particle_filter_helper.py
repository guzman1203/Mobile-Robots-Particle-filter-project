import heapq

# Adjacency list wall configuration
CELL_GRAPH = {
     1: [ 2,  6],       2: [ 1,  3],   3: [ 2,  4],   4: [ 3,  5],   5: [ 4, 10],
     6: [ 1, 11],       7: [ 8    ],   8: [ 7,  9],   9: [ 8, 10],  10: [ 5,  9],
    11: [ 6, 12, 16],  12: [11, 13],  13: [12, 14],  14: [13, 15],  15: [14, 20],
    16: [11, 17, 21],  17: [16, 18],  18: [17, 19],  19: [18, 20],  20: [15, 19, 25],
    21: [16, 22],      22: [21, 23],  23: [22, 24],  24: [23, 25],  25: [20, 24]
}

GRID_SIZE = 5

def get_environment_coords(cell_number):
    '''
    Assumes that the grid size is odd and center is 0,0
    
    Args:
        grid_size (int) : size of the gride or the number of cells it is wide and long
        cell_number (int) : the cell number for coordinates calculation
    
    Returns:
        y and x coordinate for the cell number given the grid size
    '''

    cell_number -= 1
    horizontal_depth = cell_number % GRID_SIZE
    vertical_depth = cell_number // GRID_SIZE
    x_coordinate = horizontal_depth - (GRID_SIZE//2)
    y_coordinate = (GRID_SIZE//2) - vertical_depth

    return y_coordinate, x_coordinate

def get_cardinal_between_neighbors(start_cell, neighbor_cell):
    cardinal = ""
    start_y, start_x = get_environment_coords(start_cell)
    neighbor_y, neighbor_x = get_environment_coords(neighbor_cell)
    y_displacement = neighbor_y - start_y
    x_displacement = neighbor_x - start_x
    
    if x_displacement == 1:
        cardinal = 'E'
    elif y_displacement == 1:
        cardinal = 'N'
    elif x_displacement == -1:
        cardinal = 'W'
    elif y_displacement == -1:
        cardinal = 'S'
    else:
        raise Exception("Error within get_cadinal_direction_between(). Cells are either not cardinal neighbors.")

    return cardinal

def get_all_shortest_paths(start):
    '''
    Using dijkstra's algorithm, returns the path strings in the from of a combination of ENWS movements for all cells.

    Args:
        start (int) : the cell number of the starting cell
    
    Returns:
        A dictionary of all shortest paths to each cell from the start cell.
    '''
    # TODO helper function to turn paths into their RLE forms. WWNNEEEESWWW == 2W2N4E1S3W. Note seems only better with larger grid sizes.

    priority_que = [(0, start, "")]
    paths = {node: (float('inf'), "") for node in CELL_GRAPH}
    paths[start] = (0,"")

    while priority_que:
        current_distance, current_node, current_path = heapq.heappop(priority_que)
        for neighbor in CELL_GRAPH[current_node]:
            distance = current_distance + 1
            if distance < paths[neighbor][0]:
                path = current_path + get_cardinal_between_neighbors(current_node, neighbor)
                paths[neighbor] = (distance, path)
                heapq.heappush(priority_que, (distance, neighbor, path))

    return paths

def get_shortest_path(start, goal):
    '''
    Uses dijkstra's algorithm to return the shortest path between start and goal cells.

    Ars:
        start (int) : starting cell number
        goal  (int) : goal cell number

    Returns:
        The shortest path between start and goal cells in cardinal movement format.
    '''
    # TODO improve via a .get for the case that goal is unreachable

    path = 1
    return get_all_shortest_paths(start)[goal][path]

def get_cell_map_coords(cell_number):
    row = cell_number // GRID_SIZE
    col = (cell_number - 1) % GRID_SIZE

    return row, col

def print_path(grid_size, start_cell, path, current_cell):
    '''
    Displays the path in the grid and the position given. Intended to represent the current location of the robot with it's path.
    '''

    # Generate the grid with cell numbers
    grid = [[col + row * grid_size + 1 for col in range(grid_size)] for row in range(grid_size)]

    # Convert start_cell and current_cell into grid coordinates
    sy, sx = get_cell_map_coords(start_cell)
    cy, cx = get_cell_map_coords(current_cell)

    # Add start cell to path_coords
    path_coords = [(sy, sx)]

    # Define movement map
    moves = {
        "E": ( 0, 1),
        "N": (-1, 0),
        "W": ( 0,-1),
        "S": ( 1, 0), 
    }

    # Generate the path coordinates
    y, x = sy, sx
    for direction in path:
        dy, dx = moves[direction]
        y += dy
        x += dx
        path_coords.append((y, x))

    # Calculate the width needed for formatting
    max_width = len(str(grid_size*grid_size))

    print("\nCell Map:")
    print("-|-".join(["-"*max_width] * grid_size))
    for row in range(grid_size):
        row_output = []
        for col in range(grid_size):
            if row == cy and col == cx:
                cell = "R"
            elif row == sy and col == sx:
                cell = "S"
            elif (row, col) in path_coords:
                cell = "~"
            else:
                cell = str(grid[row][col])
            row_output.append(cell.center(max_width))
        print(" | ".join(row_output) )
        print("-|-".join(["-"*max_width] * grid_size))

if __name__ == "__main__":
    print_path(5, 13, "WWNNEEEESWW", 7)