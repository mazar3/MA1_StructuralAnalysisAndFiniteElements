import numpy as np
from .elements import ElementType, tetra4_K
from provided_code.dofpos import compute_dofpos

def assembly(node: np.ndarray,
             elem: np.ndarray,
             eltp,
             mater: np.ndarray) -> np.ndarray:
    """
    [Ksys] = assembly(node, elem, eltp, mater)

    Function that assembles the structural stiffness matrix from the element
    stiffness matrix contributions.

    Inputs:
      node  - nodal coordinates matrix of the structure [nnode x 3]
      elem  - element definition matrix [nelem x (2+nen)],
              columns = [type#, material#, node1, node2, ...] (node# 1-based)
      eltp  - element type table (dict): maps type# -> type name, e.g. {1: 'tetra4'}
      mater - material properties matrix of the structure [nmat x p]
              (rows addressed by material# 1-based)

    Output:
      Ksys  - stiffness matrix of the complete system [ndof x ndof],
              where ndof = 3 * nnode

    Notes:
      - Element node IDs and material IDs are expected 1-based (MATLAB style);
        conversion to 0-based is handled internally.
      - Currently implemented element type: 'tetra4'.
    """

    # TODO: Number of elements in the structure
    nelem = 

    # TODO: Total number of degrees of freedom of the structure
    ndof  = 

    # TODO: Initialize Ksys full of zeros
    # For efficient storage, one should use a sparse matrix from scipy.sparse instead of NumPy (not necessary in this project)
    Ksys = 

    for e in range(nelem):
        type_id = int(elem[e, 0])

        # Checks if the provided type is in the enum from elements.py (no need to change this!)
        try:
            eltpe   = ElementType(eltp[type_id])
        except:
            raise NotImplementedError(f"Element type '{eltp[type_id]}' not implemented.")

        match eltpe:
            case ElementType.TETRA4:
                # You will need dofpos = compute_dofpos(nnode)

                # TODO: Determine the position of the dof in the unknowns

                # TODO: Add element contribution to the structural stiffness matrix

                # !!! Python is a 0-based language !!!
            case _:
                raise NotImplementedError(f"Element type '{eltp[type_id]}' not implemented.")


    return Ksys
