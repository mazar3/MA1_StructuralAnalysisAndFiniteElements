import numpy as np
from typing import Tuple
from provided_code.dofpos import compute_dofpos

def boundcond(pdof: np.ndarray,
              Ksys: np.ndarray,
              Fext: np.ndarray,
              method: float,
              nnode: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply Dirichlet boundary conditions to the structural system.

    Parameters
    ----------
    pdof : (npdof, 3) array
        Prescribed DOFs: [node#, dir, value] with 1-based node# and dir ∈ {1,2,3}.
    Ksys : (ndof, ndof) array
        Global stiffness matrix before BCs (dense).
    Fext : (ndof,) array
        External force vector before BCs.
    method : float
        If method < 0  → direct elimination method.
        If method >= 0 → penalty method with Z = method.
    nnode : int
        Total number of nodes in the structure.

    Returns
    -------
    Kcl, Fcl : tuple of np.ndarray
        - Kcl : Stiffness matrix with BCs applied (same shape as Ksys)
        - Fcl : Force vector with BCs applied (same shape as Fext)

    A tuple is the standard way to return multiple values in Python.
    Python functions always return a single object, and when you write
    `return Kcl, Fcl`, those values are automatically bundled into a tuple.
    """

    # TODO: Total number of degrees of freedom of the system
    ndof = 

    # Build dofpos: shape (nnode, 3)
    dofpos = compute_dofpos(nnode)

    # TODO: Number of nodal degrees of freedom
    ndof_node = 

    # TODO: Initialize Kcl and Fcl 
    # !!! Numpy arrays are mutable, so can be modified outside of the function scope. See .copy()
    Fcl = 
    Kcl = 

    if method < 0:
        # TODO: Direct method
    else:
        # TODO: Penalty method
        Z = float(method)

    return Kcl, Fcl
