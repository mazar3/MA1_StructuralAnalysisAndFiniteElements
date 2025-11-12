import numpy as np

import code_to_be_implemented.elements as elements
from code_to_be_implemented.elements import tetra4_coord_transform, tetra4_strain_stress

n = np.array([[0,0,0],
                  [1,0,0],
                  [0,1,0],
                  [0,0,1]])
m = ([1000,.3])
u = np.array([0,0,0,1,0,0,0,0,0,0,0,0])

strain, stress = tetra4_strain_stress(n, m, u)
print(strain)