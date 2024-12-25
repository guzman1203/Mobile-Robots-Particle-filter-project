Python program to compute the robot's shortest path to a Goal Cell, using Dijkstra's algorithm. 
This algorithm combines a predefined cell wall configuration and sensor inputs, such as camera data, IMU, or LiDAR. 
The program employs a data structure to model the maze, detailing cell connectivity and obstacles. 
By leveraging the predefined wall configuration, the robot identifies and avoids blocked cells on its path to the goal. 
The testing scenario involves a 5x5 grid maze, Maze 8, as shown in Figure 1. Maze 8 consists of 25 cells with wall obstacles. 
The program predefines the shortest path and total distance to the goal before the robot begins its movement. 
As the robot progresses through the maze, the program outputs the cell coordinates, cell number, and the robot’s orientation at each step.
