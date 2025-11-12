import numpy as np
from code_to_be_implemented.elements import tetra4_strain_stress

nodee = np.array([[0, 0, 0],
                  [1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]])

m = np.array([1000, 0.3])

u_cases = [
    np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),     # node2 on X
    np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]),     # node2 on Y
    np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]),     # node2 on Z
    np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]),     # node3 on X
    np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]),     # node3 on Y
    np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]),     # node3 on Z
    np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]),     # node4 on X
    np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]),     # node4 on Y
    np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])      # node4 on Z
]

good_answers = [
    np.array([1., 0., 0., 0., 0., 0.]),
    np.array([0., 0., 0., 1., 0., 0.]),
    np.array([0., 0., 0., 0., 0., 1.]),
    np.array([0., 0., 0., 1., 0., 0.]),
    np.array([0., 1., 0., 0., 0., 0.]),
    np.array([0., 0., 0., 0., 1., 0.]),
    np.array([0., 0., 0., 0., 0., 1.]),
    np.array([0., 0., 0., 0., 1., 0.]),
    np.array([0., 0., 1., 0., 0., 0.])
]


for i in range(9):
    u = u_cases[i]
    good_answer = good_answers[i]

    strain, stress = tetra4_strain_stress(nodee, m, u)

    print(f"Test {i+1}:")
    print(f"  > {strain}")
    print(f"  > {good_answer}")
    if np.allclose(strain,good_answer): print(f"    > Correct !!!!!")
    else:print(f"    > ERROR !!!!!")