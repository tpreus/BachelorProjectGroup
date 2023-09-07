# BachelorProjectGroup
This program was developed during the bachelor project "Elastic shape matching" at the University of Bonn.

The program requires:
gptoolbox (https://github.com/alecjacobson/gptoolbox)
YALMIP (https://yalmip.github.io/), download and add to matlab path
Gurobi (https://www.gurobi.com/solutions/gurobi-optimizer/), download, obtain licens and add to matlab path
Mosek (https://www.mosek.com/) download, obtain licens and add to matlab path (optional)
the feature calculation is adopted from Emanuele Rodolà

The program presents an approach to the three dimensional shape matching problem by placing points of mesh Y on the surface of mesh X. 
A mixed integer linear optimization problem with different constraints is implemented

The frame.m file is an example application to the problem

The getConstraints.m file takes the X and Y meshes and a vector indicating which constraints to compute, and returns the matrix A and vector b for optimization
The getObjectiveFunction.m file takes the meshes X and Y and their feature vectors and returns the objective vector c for the optimization.

The problem min c^t x s.t. Ax \leq b is then optimized in frame.m

getCoordinated.m takes mesh X and Y and the values x from the optimization program and returns a list of coordinates [x,y,z] for the points of Y on the surface from x

