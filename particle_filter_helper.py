import heapq
import re
from colorama import Fore, Style
from constants import GRID_SIZE, CELL_SIZE, CELL_GRAPH


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

    path = 1 # index of the path
    return get_all_shortest_paths(start)[goal][path]

def get_cell_map_coords(cell_number):
    row = cell_number // GRID_SIZE
    col = (cell_number - 1) % GRID_SIZE

    return row, col

def get_cell_number(row, column):
    return row * GRID_SIZE + column + 1

# Helper function to calculate visible length of a string
def visible_length(text):
    color_code_pattern = re.compile(r'\x1b\[[0-9;]*m')  # Matches ANSI color codes
    return len(color_code_pattern.sub('', text))  # Remove color codes and calculate length

# Helper function to center text with color codes
def color_center(text, width):
    text_length = visible_length(text)
    padding = (width - text_length) // 2
    return " " * padding + text + " " * (width - text_length - padding)
def print_path(grid_size, start_cell, path, current_cell):
    '''
    Displays the path in the grid and the position given. Intended to represent the current location of the robot with it's path.
    '''

    # Generate the grid with cell numbers
    grid = [[col + row * grid_size + 1 for col in range(grid_size)] for row in range(grid_size)]

    # Convert start_cell and current_cell into grid coordinates
    sy, sx = get_cell_map_coords(start_cell)
    cy, cx = get_cell_map_coords(current_cell)

    # Add start cell to path coordinates
    path_coords = [(sy, sx)]

    # Add start cell to path cell numbers
    path_cell_numbers = [f"{Fore.GREEN}{str(start_cell)}"]

    # Define movement map
    moves = {
        "E": ( 0, 1),
        "N": (-1, 0),
        "W": ( 0,-1),
        "S": ( 1, 0), 
    }

    # Generate the path coordinates
    row, col = sy, sx
    for direction in path:
        dy, dx = moves[direction]
        row += dy
        col += dx
        path_coords.append((row, col))
        path_cell_numbers.append(f"{Fore.BLUE}{str(get_cell_number(row, col))}")
    path_cell_numbers[-1] = f"{Fore.RED}{get_cell_number(row, col)}"


    # Calculate the width needed for formatting
    max_width = len(str(grid_size*grid_size))

    # Print cell numbers on path
    print("\nCells on path:")
    print(", ".join(path_cell_numbers))

    # Print path on cell map and current position
    print("\nCell Map:")
    print("-|-".join(["-"*max_width] * grid_size))
    for row in range(grid_size):
        row_output = []
        for col in range(grid_size):
            if row == cy and col == cx:
                cell = f"{Fore.RED}{get_cell_number(row, col)}{Style.RESET_ALL}"
            elif row == sy and col == sx:
                cell = f"{Fore.GREEN}{start_cell}{Style.RESET_ALL}"
            elif (row, col) in path_coords:
                cell = f"{Fore.BLUE}{get_cell_number(row, col)}{Style.RESET_ALL}"
            else:
                cell = str(grid[row][col])
            row_output.append(color_center(cell, max_width))
        print(" | ".join(row_output))
        print("-|-".join(["-"*max_width] * grid_size))
    print(f"{Fore.BLUE}path{Style.RESET_ALL}, {Fore.GREEN}starting cell{Style.RESET_ALL}, {Fore.RED}current location{Style.RESET_ALL}")


if __name__ == "__main__":
    # testing
    print_path(5, 13, "WWNNEEEESWWW", 7)