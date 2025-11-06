import importlib.util
import os
from types import ModuleType
import numpy as np
from .dofpos import compute_dofpos
from code_to_be_implemented.assembly import assembly
from code_to_be_implemented.boundcond import boundcond
from code_to_be_implemented.elements import ElementType, tetra4_strain_stress
from typing import Dict, Any

def linel_fem_solver(input_path: str):
    # Loads the variables from the input file to the inp module (example: the node variable from the input file is now inp.node) and validates them
    inp = load_input_module(input_path)
    validate_input_module(inp)

    fname = os.path.basename(input_path)
    print(f"\n🟢 Input '{fname}' is properly formatted.\n")

    # Show some info on the data
    nnode = inp.node.shape[0]
    nelem = inp.elem.shape[0]
    ndof  = 3 * nnode
    npdof = 0 if inp.pdof.size == 0 else inp.pdof.shape[0]
    nnodf = 0 if inp.nodf.size == 0 else inp.nodf.shape[0]

    print(f"  datafile = {fname}\n")
    print(f"  nelem = {nelem:9d}  nnode = {nnode:9d}  ndof = {ndof:9d}  npdof = {npdof:9d}  nnodf = {nnodf:9d}\n")

    # ----------------------- EXTERNAL FORCE VECTOR -----------------------
    # the external point load force vector
    Fp = getpforce(inp.nodf, nnode)

    # no volume or surface force is considered!

    Fext = Fp

    # -------------------------- STIFFNESS MATRIX -------------------------
    # assembly of the structural stiffness matrix
    Ksys = assembly(inp.node, inp.elem, inp.eltp, inp.mater)

    # apply the boundary conditions
    Kcl, Fcl = boundcond(inp.pdof, Ksys, Fext, inp.bc_method, nnode)

    # -------------------------- SOLVE THE SYSTEM -------------------------
    try:
        u = np.linalg.solve(Kcl, Fcl)
    except np.linalg.LinAlgError as e:
        raise RuntimeError(f"❌ Failed to solve linear system: {e}")
    print(f"\n✅ Input '{fname}' is solved.\n")

    # ------------------------- SAVES OUTPUT DATA -------------------------

    R = getreaction(Ksys, u, inp.node.shape[0])
    strain, stress = getstrainstress(inp.node, inp.elem, inp.eltp, inp.mater, u)
    paths = save_results(input_path, u=u.reshape(-1,3), R=R, strain=strain, stress=stress)
    print("Saved:", paths["npz"])
    print("Saved:", paths["txt"])

# ======================== Main FEM logic functions =======================

def getpforce(nodf: np.ndarray, nnode: int) -> np.ndarray:

    """
    Create the global force vector from the specified nodal point loads.

    Parameters
    ----------
    nodf : np.ndarray
        Nodal force array from input file, shape (nnodf, 3):
        [node#, dir, value], 1-based indexing for node# and dir.
    dofpos : np.ndarray
        Positioning of the degrees of freedom in the unknowns, shape (nnode, 3).
        dofpos[i, j] = global dof index (0-based) for node i, dir j.

    Returns
    -------
    Fp : np.ndarray
        Force vector corresponding to the point loads, shape (ndof,)
    """

    ndof = 3 * nnode
    # dofpos[i, j] gives global dof index (0-based) for node i and dir j
    dofpos = compute_dofpos(nnode)

    Fp = np.zeros(ndof, dtype=float)
    if nodf.size == 0:
        return Fp

    # convert 1-based node and dir to 0-based indices
    nodes = nodf[:, 0].astype(int) - 1
    dirs  = nodf[:, 1].astype(int) - 1
    vals  = nodf[:, 2].astype(float)

    # map (node, dir) -> global dof index
    gdofs = dofpos[nodes, dirs]

    # accumulate loads in case multiple loads hit the same dof
    np.add.at(Fp, gdofs, vals)
    return Fp

def getreaction(K: np.ndarray, u: np.ndarray, nnode: int) -> np.ndarray:
    """
    Compute the nodal reaction forces of the structure.

    Parameters
    ----------
    K : np.ndarray, shape (ndof, ndof)
        Structural stiffness matrix.
    u : np.ndarray, shape (ndof,) or (ndof, 1)
        Global displacement vector.
    nnode : int
        Number of nodes in the structure.

    Returns
    -------
    R : np.ndarray, shape (nnode, 3)
        Nodal reaction forces [Rx, Ry, Rz] for each node.
    """
    # Compute global reaction force vector: Freac = K * u
    Freac = K @ u.reshape(-1)  # ensure 1D

    # Build DOF-position mapping
    dofpos = compute_dofpos(nnode)

    # Rearrange into nodal format
    R = np.zeros((nnode, 3), dtype=float)
    R[:, 0] = Freac[dofpos[:, 0]]  # X-direction reactions
    R[:, 1] = Freac[dofpos[:, 1]]  # Y-direction reactions
    R[:, 2] = Freac[dofpos[:, 2]]  # Z-direction reactions

    return R

def getstrainstress(node: np.ndarray,
                    elem: np.ndarray,
                    eltp: dict,
                    mater: np.ndarray,
                    u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute per-element strain and stress.

    Parameters
    ----------
    node : (nnode, 3) array
        Nodal coordinates.
    elem : (nelem, >=6) array
        Element table: [type#, material#, n1, n2, n3, n4, ...] with 1-based IDs.
    eltp : dict
        Maps type# -> element name (e.g., {1: "tetra4"}).
    mater : (nmat, p) array
        Material properties; row (material# - 1) is the material for an element.
    u : (ndof,) or (ndof,1) array
        Global displacement vector.

    Returns
    -------
    strain : (nelem, 6) array
    stress : (nelem, 6) array
    """
    nnode = node.shape[0]
    ndof = 3 * nnode
    u_flat = u.reshape(-1)
    if u_flat.size != ndof:
        raise ValueError(f"u has length {u_flat.size}, expected {ndof} (3*nnode).")

    # dofpos[i,:] = [ux_i, uy_i, uz_i] global DOF indices
    dofpos = np.column_stack((
        np.arange(0, ndof, 3),
        np.arange(1, ndof, 3),
        np.arange(2, ndof, 3),
    ))

    nelem = elem.shape[0]
    strain = np.zeros((nelem, 6), dtype=float)
    stress = np.zeros((nelem, 6), dtype=float)

    for ie in range(nelem):
        type_id = int(elem[ie, 0])
        mat_id  = int(elem[ie, 1]) - 1  # 0-based
        etype   = eltp[type_id]

        # Checks if the provided type is in the enum from elements.py (no need to change this!)
        try:
            eltpe   = ElementType(eltp[type_id])
        except:
            raise NotImplementedError(f"Element type '{eltp[type_id]}' not implemented.")

        match eltpe:
            case ElementType.TETRA4:
                # element node ids (0-based)
                node_ids = elem[ie, 2:6].astype(int) - 1
                nodee = node[node_ids, :]                # (4,3)

                # element displacement vector ue (12,), order [u1x,u1y,u1z,u2x,...,u4z]
                ue_idx = dofpos[node_ids, :].reshape(-1) # (12,)
                ue = u_flat[ue_idx]

                # per-element strain & stress
                e6, s6 = tetra4_strain_stress(nodee, mater[mat_id, :], ue)
                strain[ie, :] = e6
                stress[ie, :] = s6
            case _:
                raise NotImplementedError(f"Element type '{eltp[type_id]}' not implemented.")

    return strain, stress

# Runs a few checks on the input, if an error occurs at here, your input file is probably the issue
def validate_input_module(inp):
    required_vars = ["node", "elem", "eltp", "mater", "pdof", "nodf", "bc_method"]

    for var in required_vars:
        if not hasattr(inp, var):
            raise AttributeError(f"Missing required variable: {var}")

    if not isinstance(inp.node, np.ndarray):
        raise TypeError(f"'node' must be a NumPy array, got {type(inp.node)}")
    if not isinstance(inp.elem, np.ndarray):
        raise TypeError(f"'elem' must be a NumPy array, got {type(inp.elem)}")
    if not isinstance(inp.mater, np.ndarray):
        raise TypeError(f"'mater' must be a NumPy array, got {type(inp.mater)}")
    if not isinstance(inp.pdof, np.ndarray):
        raise TypeError(f"'pdof' must be a NumPy array, got {type(inp.pdof)}")
    if not isinstance(inp.nodf, np.ndarray):
        raise TypeError(f"'nodf' must be a NumPy array, got {type(inp.nodf)}")
    if not isinstance(inp.eltp, dict):
        raise TypeError(f"'eltp' must be a dict, got {type(inp.eltp)}")
    if not isinstance(inp.bc_method, int):
        raise TypeError(f"'bc_method' must be an int, got {type(inp.bc_method)}")

    # shape checks
    if inp.node.ndim != 2 or inp.node.shape[1] != 3:
        raise ValueError(f"'node' must be (nnode, 3), got {inp.node.shape}")
    if inp.elem.ndim != 2 or inp.elem.shape[1] < 6:
        raise ValueError(f"'elem' must have at least 6 columns, got {inp.elem.shape}")
    if inp.pdof.size and inp.pdof.shape[1] != 3:
        raise ValueError(f"'pdof' must have 3 columns, got {inp.pdof.shape}")
    if inp.nodf.size and inp.nodf.shape[1] != 3:
        raise ValueError(f"'nodf' must have 3 columns, got {inp.nodf.shape}")

    # 1-based indexing checks
    if not np.isin(inp.elem[:, 0], list(inp.eltp.keys())).all():
        raise ValueError("'elem' type not in eltp dictionary")
    if inp.elem[:, 1].min() < 1 or inp.elem[:, 1].max() > inp.mater.shape[0]:
        raise ValueError("'elem' material not in material matrix")
    if inp.elem[:, 2:].min() < 1 or inp.elem[:, 2:].max() > inp.node.shape[0]:
        raise ValueError("'elem' node not in node matrix")
    
    if inp.pdof.size and (inp.pdof[:, 0].min() < 1 or inp.pdof[:, 0].max() > inp.node.shape[0]):
        raise ValueError("'pdof' node not in node matrix")
    if inp.pdof.size and (inp.pdof[:, 1].min() < 1 or inp.pdof[:, 1].max() > 3):
        raise ValueError("'pdof' direction must be 1, 2, or 3")

    if inp.nodf.size and (inp.nodf[:, 0].min() < 1 or inp.nodf[:, 0].max() > inp.node.shape[0]):
        raise ValueError("'nodf' node not in node matrix")
    if inp.nodf.size and (inp.nodf[:, 1].min() < 1 or inp.nodf[:, 1].max() > 3):
        raise ValueError("'nodf' direction must be 1, 2, or 3")


# ======================== DATA HANDLING =======================
# The rest of the code (below) has no interrest for the FEM course, it is just to handle data stream

# Load a .py input file and return the module so variables can be read.
def load_input_module(path: str) -> ModuleType:
    abspath = os.path.abspath(path)
    if not os.path.exists(abspath):
        raise FileNotFoundError(f"Input file not found: {abspath}")

    spec = importlib.util.spec_from_file_location("fem_input", abspath)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load Python file: {abspath}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    return mod

def save_results(input_path: str, **arrays: Any) -> Dict[str, str]:
    """
    Save arbitrary named arrays next to the input .py file:
      - <stem>_results.npz  (compressed, easy to load)
      - <stem>_results.txt  (human-readable)

    Usage:
      save_results(input_path, u=u, R=R, strain=strain, stress=stress)

    Returns
    -------
    paths : dict
        {"npz": "<..._results.npz>", "txt": "<..._results.txt>"}
    """
    if not arrays:
        raise ValueError("No data provided to save_results(). Pass named arrays, e.g. u=u.")

    in_dir  = os.path.dirname(os.path.abspath(input_path))
    stem    = os.path.splitext(os.path.basename(input_path))[0]

    # 1) Save NPZ (convert everything to np.ndarray)
    npz_path = os.path.join(in_dir, f"{stem}_results.npz")
    np.savez_compressed(npz_path, **{k: np.asarray(v) for k, v in arrays.items()})

    # 2) Save TXT (pretty print 1D/2D with savetxt; fallback for higher-D)
    txt_path = os.path.join(in_dir, f"{stem}_results.txt")
    with open(txt_path, "w") as f:
        f.write(f"# Results for {stem}\n\n")
        for name, arr in arrays.items():
            a = np.asarray(arr)
            f.write(f"## {name}  shape={a.shape}  dtype={a.dtype}\n")

            if a.ndim in (1, 2):
                # Use scientific notation, 3 decimal digits, align columns
                np.savetxt(f, a, fmt="%12.3e")
            else:
                # For 3D+ arrays, use array2string with similar formatting
                f.write(np.array2string(
                    a,
                    formatter={"float_kind": lambda x: f"{x:12.3e}"},
                    max_line_width=120
                ))
            f.write("\n\n")

    return {"npz": npz_path, "txt": txt_path}
