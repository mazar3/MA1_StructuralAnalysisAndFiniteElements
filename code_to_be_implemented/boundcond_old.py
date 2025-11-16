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
    """

    # Total number of degrees of freedom of the system
    # (we take it directly from the size of Ksys)
    ndof = Ksys.shape[0]

    # Build dofpos: shape (nnode, 3)
    # dofpos[node_idx, local_dir] = global dof index (0-based)
    dofpos = compute_dofpos(nnode)

    # Number of nodal degrees of freedom
    ndof_node = dofpos.shape[1]   # normally 3 in 3D

    # Initialize Kcl and Fcl (work on copies to avoid modifying inputs)
    Fcl = Fext.copy()
    Kcl = Ksys.copy()

    # Number of prescribed DOFs
    npdof = pdof.shape[0]

    if method < 0:
        # ============================
        # Direct elimination method
        # ============================
        for i in range(npdof):
            # Extract from pdof: node#, dir, value
            node_id = int(pdof[i, 0])   # 1-based
            loc_dir = int(pdof[i, 1])   # 1-based (1=x, 2=y, 3=z)
            value   = float(pdof[i, 2])

            # Convert to 0-based indices for Python
            node_idx = node_id - 1
            dir_idx  = loc_dir - 1

            # Global degree of freedom index
            dof = int(dofpos[node_idx, dir_idx])

            # Apply direct method:
            # We enforce q[dof] = value
            # K * q = F  ->  modify K and F so that the solution satisfies this
            Fcl = Fcl - Kcl[:, dof] * value   # move known term to RHS

            # Zero out the column and row, then put 1 on the diagonal
            Kcl[:, dof] = 0.0
            Kcl[dof, :] = 0.0
            Kcl[dof, dof] = 1.0

            # Set corresponding RHS entry to the prescribed value
            Fcl[dof] = value

    else:
        # ============================
        # Penalty method
        # ============================
        Z = float(method)   # penalty parameter

        for i in range(npdof):
            node_id = int(pdof[i, 0])   # 1-based
            loc_dir = int(pdof[i, 1])   # 1-based
            value   = float(pdof[i, 2])

            node_idx = node_id - 1
            dir_idx  = loc_dir - 1

            dof = int(dofpos[node_idx, dir_idx])

            # Add huge stiffness on the diagonal and corresponding RHS
            Kcl[dof, dof] += Z
            Fcl[dof]      += Z * value

    return Kcl, Fcl
