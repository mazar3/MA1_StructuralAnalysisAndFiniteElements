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

    # get the total number of degrees of freedom of the system
    ndof = 3 * nnode

    # build dofpos
    dofpos = compute_dofpos(nnode)

    # get the number of nodal degrees of freedom
    ndof_node = 3

    # initialize Kcl and Fcl (Numpy arrays are mutable, so can be modified outside the function scope)
    Fcl = Fext.copy()
    Kcl = Ksys.copy()

    if method < 0:  # DIRECT METHOD
        # STEP 1 (cf slides): Group known terms on RHS
        for i in range(pdof.shape[0]):
            # get node number, direction and value
            node_num = int(pdof[i,0]) - 1  # 0-based index
            node_dir = int(pdof[i,1]) - 1  # 0-based index
            node_val = pdof[i,2]

            # get the global dof position in the matrix
            node_pos = dofpos[node_num,node_dir]

            # group known terms on RHS
            Fcl = Fcl - Ksys[:,node_pos] * node_val

        # STEP 2 (cf slides): "Remove" equation corresponding to the known quantity
        for i in range(pdof.shape[0]):
            # get node number, direction and value
            node_num = int(pdof[i,0]) - 1  # 0-based index
            node_dir = int(pdof[i,1]) - 1  # 0-based index
            node_val = pdof[i,2]

            # get the global dof position in the matrix
            node_pos = dofpos[node_num, node_dir]

            # edit Kcl (put 0s to the row and column)
            Kcl[node_pos,:] = 0
            Kcl[:,node_pos] = 0
            # put 1 to the diagonal element of Kcl
            Kcl[node_pos,node_pos] = 1
            # put the value to the element of Fcl
            Fcl[node_pos] = node_val

    else: # Penalty method
        Z = method
        # loop over each prescribed degree of freedom
        for i in range(pdof.shape[0]):
            # get node number, direction, and value
            node_num = int(pdof[i,0]) - 1 # 0-based index
            node_dir = int(pdof[i,1]) - 1 # 0-based index
            node_val = pdof[i, 2]

            # get the global dof position in the matrix
            node_pos = dofpos[node_num, node_dir]

            # add the penalty term to the diagonal of the stiffness matrix
            Kcl[node_pos, node_pos] += Z

            # add the penalty term to the force vector
            Fcl[node_pos] = Fcl[node_pos] + Z*node_val

    return Kcl, Fcl