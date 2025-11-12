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

    C_t = np.vstack((np.ones(4),nodee.T))

    return C_t

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
    Ve = np.linalg.det(Ct) / 6.0
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

    Ct = tetra4_coord_transform(nodee)
    if np.isclose(np.linalg.det(Ct), 0):
        raise np.linalg.LinAlgError("Element with zero or negative volume")

    Ct_inv = np.linalg.inv(Ct)

    B = np.zeros((6, 12))

    for i in range(4):
        col = i * 3

        dxi_dx = Ct_inv[i, 1]
        dxi_dy = Ct_inv[i, 2]
        dxi_dz = Ct_inv[i, 3]

        B[0, col] = dxi_dx
        B[1, col + 1] = dxi_dy
        B[2, col + 2] = dxi_dz

        B[3, col] = dxi_dy
        B[3, col + 1] = dxi_dx

        B[4, col + 1] = dxi_dz
        B[4, col + 2] = dxi_dy

        B[5, col] = dxi_dz
        B[5, col + 2] = dxi_dx

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

    C1 = E * (1 - v) / ((1 + v) * (1 - 2 * v))
    C2 = E * v / ((1 + v) * (1 - 2 * v))
    C3 = E / (2 * (1 + v))

    H[0, 0] = C1
    H[1, 1] = C1
    H[2, 2] = C1

    H[0, 1] = C2
    H[0, 2] = C2
    H[1, 0] = C2
    H[1, 2] = C2
    H[2, 0] = C2
    H[2, 1] = C2

    H[3, 3] = C3
    H[4, 4] = C3
    H[5, 5] = C3

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

    Ke = B.T @ H @ B * Ve
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

    H = hooke_iso_3d(E, v)
    B = tetra4_B(nodee)

    strain = B @ ue
    stress = H @ strain

    return strain, stress
