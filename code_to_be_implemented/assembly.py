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
    # get the number of elements in the structure:
    nelem = elem.shape[0]

    # get the number of degrees of freedom of the structure:
    ndof  = node.shape[0]*3 #because there are 3ndof per node (x,y,z)

    # get the number of nodes:
    nnode = node.shape[0]

    # get the position of the dof in the unknowns
    dofpos = compute_dofpos(nnode)

    Ksys = np.zeros((ndof,ndof)) # For efficient storage, one should use a sparse matrix from scipy.sparse instead of NumPy (not necessary in this project)

    for e in range(nelem):
        type_id = int(elem[e, 0])

        # Checks if the provided type is in the enum from elements.py (no need to change this!)
        try:
            eltpe   = ElementType(eltp[type_id])
        except:
            raise NotImplementedError(f"Element type '{eltp[type_id]}' not implemented.")

        match eltpe:
            case ElementType.TETRA4:
                # get the properties of the element
                e_type = elem[e,1] - 1 # 0-based index
                e_material = mater[e_type,:]

                # get the 4 nodes of the element e and their position
                node_ids = elem[e,2:6] - 1 # 0-based index
                node_pos = node[node_ids,:]

                # get the stiffness matrix (12x12) of the element
                Ke = tetra4_K(node_pos, e_material)

                # Add element contribution to the structural stiffness matrix
                for row in range(12):
                    # get the local node (0,1,2 or 3) and the direction (0,1 or 2)
                    n_local_row = row // 3
                    dir_row = row % 3
                    # get the global node
                    n_global_row = node_ids[n_local_row]
                    # get the global row index in Ksys
                    new_row = dofpos[n_global_row, dir_row]

                    for col in range(12):
                        # same thing but for the columns
                        n_local_col = col // 3
                        dir_col = col % 3
                        n_global_col = node_ids[n_local_col]
                        new_col = dofpos[n_global_col, dir_col]

                        # add the contribution of the element e to Ksys
                        Ksys[new_row, new_col] = Ksys[new_row, new_col] + Ke[row, col]
            case _:
                raise NotImplementedError(f"Element type '{eltp[type_id]}' not implemented.")

    return Ksys