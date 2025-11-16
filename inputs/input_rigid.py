import numpy as np

# nodal coordinates [x, y, z]
node = np.array([
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 1, 1],
], dtype=float)

# mesh connectivity [type#, material#, node1 node2 node3 node4]
elem = np.array([
    [1, 1, 1, 2, 3, 4],
    [1, 1, 5, 2, 4, 3],
], dtype=int)

# prescribed degrees of freedom (displacements) [node# dir value]
pdof = np.array([
    [1, 1, 0.5],
    [1, 2, 0.1],
    [1, 3, 0.2],
    [2, 1, 0.5],
    [2, 2, 0.1],
    [2, 3, 0.2],
    [3, 1, 0.5],
    [3, 2, 0.1],
    [3, 3, 0.2],
    [4, 1, 0.5],
    [4, 2, 0.1],
    [4, 3, 0.2],
    [5, 1, 0.5],
    [5, 2, 0.1],
    [5, 3, 0.2],
], dtype=float)

# applied nodal loads [node# dir value]
nodf = np.array([], dtype=float)

# finite element type table
eltp = {1: "tetra4"}

# material parameters
E  = 1000.0
nu = 0.3

# material parameter matrix (rows = materials)
mater = np.array([
    [E, nu],
], dtype=float)

# Boundary condition method
# < 0 direct, >= 0 penalty
bc_method = -1