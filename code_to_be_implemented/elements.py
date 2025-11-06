from enum import Enum
import numpy as np

class ElementType(Enum):
    TETRA4 = "tetra4"

# ============================================================
# Core helpers
# ============================================================
# I created these for ease of implementation, but you are not requiered to use helper functions


# def tetra4_coord_transform(nodee: np.ndarray) -> np.ndarray:
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


# def tetra4_volume(nodee: np.ndarray) -> float:
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


# def tetra4_B(nodee: np.ndarray) -> np.ndarray:
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


# def hooke_iso_3d(E: float, v: float) -> np.ndarray:
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
    # TODO: Compute the stiffness matrix of the finite element
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
    # TODO: Post-processing step, compute the strain and stress matrix
    return strain, stress
