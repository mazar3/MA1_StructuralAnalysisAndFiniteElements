# provided_code/plotting.py
import os
import importlib.util
import numpy as np

# --- PyVista setup ---
import pyvista as pv


# ---------- path + loading helpers ----------

def _case_paths(case_name: str):
    """
    Resolve paths from a case name (e.g., 'example1' or 'example1.py').
    Returns (out_dir, stem, input_py, results_npz).
    """
    stem = os.path.splitext(os.path.basename(case_name))[0]
    out_dir = os.path.join("outputs", stem)
    input_py = os.path.join(out_dir, f"{stem}.py")
    results_npz = os.path.join(out_dir, f"{stem}_results.npz")
    return out_dir, stem, input_py, results_npz


def _load_input_module(input_py: str):
    """Import the input .py file as a module (contains node, elem, eltp, mater, ...)."""
    if not os.path.exists(input_py):
        raise FileNotFoundError(f"Input file not found: {input_py}")
    spec = importlib.util.spec_from_file_location("plot_input_module", input_py)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    return mod


def _load_results_npz(results_npz: str):
    """Load the saved results npz: keys may include u, R, strain, stress."""
    if not os.path.exists(results_npz):
        raise FileNotFoundError(f"Results file not found: {results_npz}")
    return np.load(results_npz)


# ---------- geometry helpers ----------

def _faces_tetra4():
    # local faces (0..3 local node ids)
    # (ordering chosen to be consistent; orientation doesn't matter for PyVista)
    return np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]], dtype=int)


def _u_to_U(u: np.ndarray) -> np.ndarray:
    """(3*nnode,) or (3*nnode,1) -> (nnode,3)."""
    u = np.asarray(u).reshape(-1)
    if u.size % 3 != 0:
        raise ValueError(f"u length {u.size} is not divisible by 3.")
    return u.reshape(-1, 3)


def _boundary_face_indices_and_owner(elem: np.ndarray):
    """
    Return boundary face node indices (nf,3) and the owning element index for each face.
    Faces are identified by sorted global node IDs; only faces seen once are kept.
    """
    faces_local = _faces_tetra4()
    face_count = {}
    face_owner = {}
    for e in range(elem.shape[0]):
        conn = elem[e, 2:6].astype(int) - 1  # 1-based -> 0-based global ids
        for f in faces_local:
            fnodes = conn[f]
            key = tuple(sorted(fnodes.tolist()))
            if key in face_count:
                face_count[key] += 1
            else:
                face_count[key] = 1
                face_owner[key] = (e, tuple(fnodes))

    # keep faces seen exactly once
    faces = []
    owners = []
    for key, cnt in face_count.items():
        if cnt == 1:
            e, orig = face_owner[key]
            faces.append(orig)    # keep the original triplet (orientation ok for rendering)
            owners.append(e)
    if len(faces) == 0:
        return np.empty((0,3), dtype=int), np.empty((0,), dtype=int)
    return np.asarray(faces, dtype=int), np.asarray(owners, dtype=int)


def _polydata_from_faces(points: np.ndarray, faces_ijk: np.ndarray) -> pv.PolyData:
    """
    Build a PyVista PolyData from points (n,3) and triangle indices (m,3).
    PyVista expects the 'faces' array flattened with counts: [3, i, j, k, 3, ...]
    """
    if faces_ijk.size == 0:
        # empty mesh is allowed; create minimal PolyData
        return pv.PolyData(points.copy())
    flat = np.hstack([np.insert(tri, 0, 3) for tri in faces_ijk])  # prepend 3 to each tri
    return pv.PolyData(points.copy(), flat.astype(np.int64))


# ---------- plotting primitives (PyVista) ----------

def _make_plotter(show: bool):
    """
    Create a Plotter. If not showing (headless / CI), enable off_screen to render without a window.
    """
    return pv.Plotter(off_screen=not show)


def plot_mesh(case_name: str, show: bool = True):
    """
    Plot the undeformed tetra4 mesh using inputs/<case>.py (no results required).
    Boundary faces only.
    """
    stem = os.path.splitext(os.path.basename(case_name))[0]
    input_py = os.path.join("inputs", f"{stem}.py")
    inp = _load_input_module(input_py)

    node, elem = np.asarray(inp.node, float), np.asarray(inp.elem, int)
    faces_ijk, _ = _boundary_face_indices_and_owner(elem)
    surf = _polydata_from_faces(node, faces_ijk)

    plotter = _make_plotter(show)
    plotter.add_mesh(surf, color="lightgray", show_edges=True)
    plotter.add_axes()
    plotter.set_background("white")
    return plotter, surf


def plot_deformed_mesh(case_name: str, scale: float = 0.0, show: bool = True):
    """
    Plot the (optionally deformed) tetra4 mesh (boundary faces only).
    """
    out_dir, stem, input_py, results_npz = _case_paths(case_name)
    inp = _load_input_module(input_py)
    node, elem = np.asarray(inp.node, float), np.asarray(inp.elem, int)

    node_plot = node.copy()
    if abs(scale) > 0:
        res = _load_results_npz(results_npz)
        if "u" not in res:
            raise KeyError(f"'u' not found in {results_npz}")
        U = _u_to_U(res["u"])
        node_plot = node + scale * U

    faces_ijk, _ = _boundary_face_indices_and_owner(elem)
    surf = _polydata_from_faces(node_plot, faces_ijk)

    plotter = _make_plotter(show)
    plotter.add_mesh(surf, color="lightgray", show_edges=True)
    plotter.add_axes()
    plotter.set_background("white")
    return plotter, surf


def plot_displacement(case_name: str, component: str = "mag", scale: float = 1.0, show: bool = True):
    """
    Plot nodal displacement field with a colorbar (optionally deformed),
    coloring boundary faces with per-vertex (nodal) scalars.
    """
    out_dir, stem, input_py, results_npz = _case_paths(case_name)
    inp = _load_input_module(input_py)
    node, elem = np.asarray(inp.node, float), np.asarray(inp.elem, int)

    res = _load_results_npz(results_npz)
    if "u" not in res:
        raise KeyError(f"'u' not found in {results_npz}")
    U = _u_to_U(res["u"])
    node_def = node + scale * U

    # nodal scalar field s (|u| or component)
    if component.lower() == "mag":
        s = np.linalg.norm(U, axis=1)
        title = f"|u| (scale={scale:g})"
    else:
        idx = {"ux": 0, "uy": 1, "uz": 2}[component.lower()]
        s = U[:, idx]
        title = f"{component} (scale={scale:g})"

    faces_ijk, _ = _boundary_face_indices_and_owner(elem)
    surf = _polydata_from_faces(node_def, faces_ijk)

    # attach point scalars for smooth per-vertex coloring
    surf.point_data.clear()
    surf.point_data["displ"] = s

    plotter = _make_plotter(show)
    plotter.add_mesh(
        surf,
        scalars="displ",
        cmap="viridis",
        show_edges=True,
        nan_color="gray",
        lighting=False,
        smooth_shading=False,  # set True if you want Gouraud shading
        show_scalar_bar=False,
    )
    plotter.add_scalar_bar(title="displacement", vertical=True)
    plotter.add_text(f"{stem} – displacement {title}", font_size=10)
    plotter.add_axes()
    plotter.set_background("white")
    return plotter, surf


def _von_mises_from_stress(s6: np.ndarray) -> np.ndarray:
    """Compute von Mises from stress Voigt (σxx,σyy,σzz,τyz,τxz,τxy) per element."""
    sx, sy, sz, tyz, txz, txy = (s6[:,0], s6[:,1], s6[:,2], s6[:,3], s6[:,4], s6[:,5])
    return np.sqrt(0.5*((sx-sy)**2 + (sy-sz)**2 + (sz-sx)**2) + 3.0*(tyz**2 + txz**2 + txy**2))


def plot_strain(case_name: str, component: str = "exx", scale: float = 1.0, show: bool = True):
    """
    Plot element-wise strain with colorbar (optionally deformed),
    coloring boundary faces by the element value (cell scalars).
    """
    out_dir, stem, input_py, results_npz = _case_paths(case_name)
    inp = _load_input_module(input_py)
    node, elem = np.asarray(inp.node, float), np.asarray(inp.elem, int)

    res = _load_results_npz(results_npz)
    if "u" not in res or "strain" not in res:
        raise KeyError(f"'u' and 'strain' must be in {results_npz}")
    U = _u_to_U(res["u"])
    node_def = node + scale * U

    comp_map = {"exx":0, "eyy":1, "ezz":2, "gyz":3, "gxz":4, "gxy":5}
    ei = comp_map[component.lower()]
    vals_elem = np.asarray(res["strain"])[:, ei]  # (nelem,)

    faces_ijk, owners = _boundary_face_indices_and_owner(elem)
    surf = _polydata_from_faces(node_def, faces_ijk)

    # cell scalars: one value per boundary face, coming from its owning element
    vals_faces = vals_elem[owners]
    surf.cell_data.clear()
    surf.cell_data["strain"] = vals_faces

    plotter = _make_plotter(show)
    plotter.add_mesh(
        surf,
        scalars="strain",
        cmap="viridis",
        show_edges=True,
        nan_color="gray",
        lighting=False,
        show_scalar_bar=False,
    )
    plotter.add_scalar_bar(title=f"strain {component}", vertical=True)
    plotter.add_text(f"{stem} – strain {component} (scale={scale:g})", font_size=10)
    plotter.add_axes()
    plotter.set_background("white")
    return plotter, surf


def plot_stress(case_name: str, component: str = "vm", scale: float = 1.0, show: bool = True):
    """
    Plot element-wise stress with colorbar (optionally deformed),
    coloring boundary faces by the element value (cell scalars).
    """
    out_dir, stem, input_py, results_npz = _case_paths(case_name)
    inp = _load_input_module(input_py)
    node, elem = np.asarray(inp.node, float), np.asarray(inp.elem, int)

    res = _load_results_npz(results_npz)
    if "u" not in res or "stress" not in res:
        raise KeyError(f"'u' and 'stress' must be in {results_npz}")
    U = _u_to_U(res["u"])
    node_def = node + scale * U

    s6 = np.asarray(res["stress"])  # (nelem,6)
    if component.lower() == "vm":
        vals_elem = _von_mises_from_stress(s6)
        label = "von Mises"
    else:
        comp_map = {"sxx":0, "syy":1, "szz":2, "tyz":3, "txz":4, "txy":5}
        vals_elem = s6[:, comp_map[component.lower()]]
        label = component

    faces_ijk, owners = _boundary_face_indices_and_owner(elem)
    surf = _polydata_from_faces(node_def, faces_ijk)

    vals_faces = vals_elem[owners]
    surf.cell_data.clear()
    surf.cell_data["stress"] = vals_faces

    plotter = _make_plotter(show)
    plotter.add_mesh(
        surf,
        scalars="stress",
        cmap="viridis",
        show_edges=True,
        nan_color="gray",
        lighting=False,
        show_scalar_bar=False,
    )
    plotter.add_scalar_bar(title=f"stress {label}", vertical=True)
    plotter.add_text(f"{stem} – stress {label} (scale={scale:g})", font_size=10)
    plotter.add_axes()
    plotter.set_background("white")
    return plotter, surf


# ---------- CLI ----------

if __name__ == "__main__":
    import argparse, sys

    parser = argparse.ArgumentParser(description="FEM plotting (PyVista)")
    parser.add_argument("case", help="case name, e.g. 'input_example'")
    parser.add_argument("--what", choices=["mesh", "displ", "strain", "stress", "defmesh"],
                        default="mesh", help="what to plot")
    parser.add_argument("--component", default="", help=(
        "for displ: one of ux,uy,uz,mag (default=mag); "
        "for strain: exx,eyy,ezz,gyz,gxz,gxy; "
        "for stress: vm,sxx,syy,szz,tyz,txz,txy"
    ))
    parser.add_argument("--scale", type=float, default=1.0,
                        help="deformation scale factor (ignored for mesh)")
    parser.add_argument("--no-show", action="store_true",
                        help="render off-screen (no window) – useful when saving figures")
    parser.add_argument("--save", default="", help="path to save .png (optional)")
    args = parser.parse_args()

    try:
        show_flag = not args.no_show

        if args.what == "mesh":
            plotter, _ = plot_mesh(args.case, show=show_flag)
        elif args.what == "displ":
            comp = args.component or "mag"
            plotter, _ = plot_displacement(args.case, component=comp, scale=args.scale, show=show_flag)
        elif args.what == "strain":
            comp = args.component or "exx"
            plotter, _ = plot_strain(args.case, component=comp, scale=args.scale, show=show_flag)
        elif args.what == "stress":
            comp = args.component or "vm"
            plotter, _ = plot_stress(args.case, component=comp, scale=args.scale, show=show_flag)
        elif args.what == "defmesh":
            plotter, _ = plot_deformed_mesh(args.case, scale=args.scale, show=show_flag)
        else:
            raise ValueError("Unknown plot type")

        # Render and optionally save
        if args.save:
            # If showing, PyVista will open a window and still write the screenshot.
            plotter.show(screenshot=args.save, auto_close=True)
            print(f"Saved: {args.save}")
        else:
            if show_flag:
                plotter.show(auto_close=True)
            else:
                # Off-screen render without saving – do a dry render
                plotter.show(auto_close=True)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
