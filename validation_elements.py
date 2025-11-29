import numpy as np
from code_to_be_implemented.elements import tetra4_strain_stress

nodee = np.array([[0, 0, 0],
                  [1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]])

m = np.array([1000, 0.3])

u_cases = [
    np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),         # Test zero
    np.array([1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0]),         # Trans X
    np.array([0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0]),         # Trans Y
    np.array([0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1]),         # Trans Z
    np.array([0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0]),        # Rot X
    np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0]),        # Rot Y
    np.array([0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0]),        # Rot Z
    np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),         # exx = 1
    np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]),         # eyy = 1
    np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]),         # ezz = 1
    np.array([0, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0, 0, 0]),     # gxy = 1
    np.array([0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.5, 0]),     # gyz = 1
    np.array([0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0])      # gxz = 1
]

good_strains = [
    np.array([0., 0., 0., 0., 0., 0.]),     # Test zero
    np.array([0., 0., 0., 0., 0., 0.]),     # Trans X
    np.array([0., 0., 0., 0., 0., 0.]),     # Trans Y
    np.array([0., 0., 0., 0., 0., 0.]),     # Trans Z
    np.array([0., 0., 0., 0., 0., 0.]),     # Rot X
    np.array([0., 0., 0., 0., 0., 0.]),     # Rot Y
    np.array([0., 0., 0., 0., 0., 0.]),     # Rot Z
    np.array([1., 0., 0., 0., 0., 0.]),     # exx
    np.array([0., 1., 0., 0., 0., 0.]),     # eyy
    np.array([0., 0., 1., 0., 0., 0.]),     # ezz
    np.array([0., 0., 0., 1., 0., 0.]),     # gxy
    np.array([0., 0., 0., 0., 1., 0.]),     # gyz
    np.array([0., 0., 0., 0., 0., 1.])      # gxz
]

E = m[0]
v = m[1]
C1 = E*(1-v) / ((1+v)*(1-2*v))
C2 = E*v / ((1+v)*(1-2*v))
C3 = E / (2*(1+v))

good_stresses = [
    np.array([0., 0., 0., 0., 0., 0.]),     # Test zero
    np.array([0., 0., 0., 0., 0., 0.]),     # Trans X
    np.array([0., 0., 0., 0., 0., 0.]),     # Trans Y
    np.array([0., 0., 0., 0., 0., 0.]),     # Trans Z
    np.array([0., 0., 0., 0., 0., 0.]),     # Rot X
    np.array([0., 0., 0., 0., 0., 0.]),     # Rot Y
    np.array([0., 0., 0., 0., 0., 0.]),     # Rot Z
    np.array([C1, C2, C2, 0., 0., 0.]),     # exx
    np.array([C2, C1, C2, 0., 0., 0.]),     # eyy
    np.array([C2, C2, C1, 0., 0., 0.]),     # ezz
    np.array([0., 0., 0., C3, 0., 0.]),     # gxy
    np.array([0., 0., 0., 0., C3, 0.]),     # gyz
    np.array([0., 0., 0., 0., 0., C3])      # gxz
]

for i in range(13):
    u = u_cases[i]
    good_strain = good_strains[i]
    good_stress = good_stresses[i]

    strain, stress = tetra4_strain_stress(nodee, m, u)

    print(f"\nTest {i + 1}:")

    print(f"  Strain: {strain}")
    if np.allclose(strain, good_strain):
        print(f"    > OK")
    else:
        print(f"    > ERROR !! (good answer : {good_strain})")

    print(f"  Stress: {stress}")
    if np.allclose(stress, good_stress):
        print(f"    > OK")
    else:
        print(f"    > ERROR !! (good answer : {good_stress})")