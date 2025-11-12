from enum import Enum
import numpy as np

class ElementType(Enum):
    TETRA4 = "tetra4"

# ============================================================
# Core helpers
# ============================================================
# I created these for ease of implementation, but you are not requiered to use helper functions


def tetra4_coord_transform(nodee: np.ndarray) -> np.ndarray:
    """
    Build the 4×4 coordinate transformation matrix for a 4-node tetrahedron.

    Parameters
    ----------
    nodee : np.ndarray, shape (4, 3)
        Element nodal coordinates. Each row is [x, y, z].

    Returns
    -------
    Ct : np.ndarray, shape (4, 4)
        Coordinate transform matrix:
            [ 1   1   1   1
              x1  x2  x3  x4
              y1  y2  y3  y4
              z1  z2  z3  z4 ]
    """

    Ct = np.zeros((4,4))
    Ct[0,:] = 1
    Ct[1:4,:] = nodee.T

    return Ct

def tetra4_volume(nodee: np.ndarray) -> float:
    """
    Compute the volume of a 4-node tetrahedron.

    Parameters
    ----------
    nodee : np.ndarray, shape (4, 3)
        Element nodal coordinates.

    Returns
    -------
    Ve : float
        Element volume.
    """

    Ct = tetra4_coord_transform(nodee)
    Ve = np.linalg.det(Ct)/6

    return Ve


def tetra4_B(nodee: np.ndarray) -> np.ndarray:
    """
    Build the 6×12 strain–displacement matrix for a 4-node tetrahedron.

    Parameters
    ----------
    nodee : np.ndarray, shape (4, 3)
        Element nodal coordinates.

    Returns
    -------
    B : np.ndarray, shape (6, 12)
        Strain–displacement matrix such that strain = B @ u_e.
    """

    B = np.zeros((6,12))
    Ve = tetra4_volume(nodee)

    x = nodee[:,0]
    y = nodee[:,1]
    z = nodee[:,2]

    a1 = (y[3]-y[1])*(z[2]-z[1]) - (y[2]-y[1])*(z[3]-z[1])
    b1 = (x[2]-x[1])*(z[3]-z[1]) - (x[3]-x[1])*(z[2]-z[1])
    c1 = (x[3]-x[1])*(y[2]-y[1]) - (x[2]-x[1])*(y[3]-y[1])

    a2 = (y[2]-y[0])*(z[3]-z[2]) - (y[3]-y[2])*(z[2]-z[0])
    b2 = (x[3]-x[2])*(z[2]-z[0]) - (x[2]-x[0])*(z[3]-z[2])
    c2 = (x[2]-x[0])*(y[3]-y[2]) - (x[3]-x[2])*(y[2]-y[0])

    a3 = (y[3]-y[0])*(z[1]-z[3]) - (y[1]-y[3])*(z[3]-z[0])
    b3 = (x[1]-x[3])*(z[3]-z[0]) - (x[3]-x[0])*(z[1]-z[3])
    c3 = (x[3]-x[0])*(y[1]-y[3]) - (x[1]-x[3])*(y[3]-y[0])

    a4 = (y[1]-y[0])*(z[2]-z[1]) - (y[2]-y[1])*(z[1]-z[0])
    b4 = (x[2]-x[1])*(z[1]-z[0]) - (x[1]-x[0])*(z[2]-z[1])
    c4 = (x[1]-x[0])*(y[2]-y[1]) - (x[2]-x[1])*(y[1]-y[0])

    B[0,0] = a1
    B[1,1] = b1
    B[2,2] = c1
    B[3,0] = b1
    B[3,1] = a1
    B[4,1] = c1
    B[4,2] = b1
    B[5,0] = c1
    B[5,2] = a1

    B[0,3] = a2
    B[1,4] = b2
    B[2,5] = c2
    B[3,3] = b2
    B[3,4] = a2
    B[4,4] = c2
    B[4,5] = b2
    B[5,3] = c2
    B[5,5] = a2

    B[0,6] = a3
    B[1,7] = b3
    B[2,8] = c3
    B[3,6] = b3
    B[3,7] = a3
    B[4,7] = c3
    B[4,8] = b3
    B[5,6] = c3
    B[5,8] = a3

    B[0,9]  = a4
    B[1,10] = b4
    B[2,11] = c4
    B[3,9]  = b4
    B[3,10] = a4
    B[4,10] = c4
    B[4,11] = b4
    B[5,9]  = c4
    B[5,11] = a4

    B = B/(6*Ve)

    return B

def hooke_iso_3d(E: float, v: float) -> np.ndarray:
    """
    Build the 6×6 Hooke matrix for 3D isotropic elasticity.

    Parameters
    ----------
    E : float
        Young's modulus.
    v : float
        Poisson's ratio.

    Returns
    -------
    H : np.ndarray, shape (6, 6)
        Constitutive matrix.
    """

    H = np.zeros((6, 6))

    C1 = E*(1-v) / ((1+v)*(1-2*v))
    C2 = E*v / ((1+v)*(1-2*v))
    C3 = E / (2*(1+v))

    H[0,0] = C1
    H[1,1] = C1
    H[2,2] = C1

    H[0,1] = C2
    H[0,2] = C2
    H[1,0] = C2
    H[1,2] = C2
    H[2,0] = C2
    H[2,1] = C2

    H[3,3] = C3
    H[4,4] = C3
    H[5,5] = C3

    return H

# ============================================================
# Public element routines
# ============================================================

def tetra4_K(nodee: np.ndarray, matere: np.ndarray) -> np.ndarray:
    """
    Compute the stiffness matrix for a 4-node tetrahedron.

    Parameters
    ----------
    nodee : np.ndarray, shape (4, 3)
        Element nodal coordinates.
    matere : np.ndarray
        Material properties [E, v].

    Returns
    -------
    Ke : np.ndarray, shape (12, 12)
        Element stiffness matrix.
    """

    E = matere[0]
    v = matere[1]

    H = hooke_iso_3d(E, v)
    B = tetra4_B(nodee)
    Ve = tetra4_volume(nodee)

    Ke = B.T @ H @ B*Ve
    return Ke


def tetra4_strain_stress(nodee: np.ndarray,
                         matere: np.ndarray,
                         ue: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute constant strain and stress for a 4-node tetrahedron.

    Parameters
    ----------
    nodee : np.ndarray, shape (4, 3)
        Element nodal coordinates.
    matere : np.ndarray
        Material properties [E, v].
    ue : np.ndarray, shape (12,)
        Element displacement vector.

    Returns
    -------
    strain : np.ndarray, shape (6,)
        Engineering strain vector.
    stress : np.ndarray, shape (6,)
        Cauchy stress vector.

    A tuple is the standard way to return multiple values in Python.
    Python functions always return a single object, and when you write
    `return Kcl, Fcl`, those values are automatically bundled into a tuple.
    """

    E = matere[0]
    v = matere[1]

    H = hooke_iso_3d(E,v)
    B = tetra4_B(nodee)

    strain = B @ ue
    stress = H @ strain

    return strain, stress
