import numpy as np

def compute_dofpos(nnode: int) -> np.ndarray:
    """
    Return the (nnode, 3) DOF-positioning array for 3D problems.
    Row i gives the global DOF indices [ux, uy, uz] for node i (0-based).
    """
    ndof = 3 * nnode
    return np.column_stack((
        np.arange(0, ndof, 3),
        np.arange(1, ndof, 3),
        np.arange(2, ndof, 3),
    ))