# ocelot/gui/layout_3d.py
#
# 3D facility plotting for MachineLayout using Plotly.
# Key improvements vs a quick prototype:
#   1) Correct metric scaling: X and Y are rendered with equal scale (round tubes stay round).
#   2) Uses W_start (entrance frame) for rigid-body placement of element meshes.
#   3) Uses explicit triangle faces for boxes/cylinders (no alphahull artifacts).
#   4) Bend “box” is aligned to the chord direction using rot_ang = -angle/2 (your fix).
#
# Expected survey item keys (from MagneticLattice.survey()):
#   item["r_start"] : np.array(3,)
#   item["r_end"]   : np.array(3,)
#   item["W_start"] : np.array(3,3)   (recommended)
#   item["W"]       : np.array(3,3)   (fallback)
#   item["element"] : element instance with .l, .width, .height, .color, .tilt, .angle

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go


# -----------------------------
# Small rotation helpers
# -----------------------------
def _rot_y(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0.0, s],
                     [0.0, 1.0, 0.0],
                     [-s, 0.0, c]], dtype=float)


def _get_roll_matrix(psi: float) -> np.ndarray:
    """
    Roll around local s-axis. Uses ocelot.common.math_op.get_tilt_matrix if available,
    otherwise falls back to a standard roll matrix around z(s)-axis in local coords
    (local axes are [x,y,s], so roll is around index 2).
    """
    try:
        from ocelot.common.math_op import get_tilt_matrix  # type: ignore
        return get_tilt_matrix(psi)
    except Exception:
        c, s = np.cos(psi), np.sin(psi)
        # roll around local s-axis (3rd axis)
        return np.array([[c, -s, 0.0],
                         [s,  c, 0.0],
                         [0.0, 0.0, 1.0]], dtype=float)


# -----------------------------
# Mesh primitives (local frame)
# -----------------------------
def _box_mesh_local(length: float, width: float, height: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Box aligned with local axes:
      x: width, y: height, s: length, with s in [0, length]
    Returns:
      verts (N,3), i, j, k arrays for Mesh3d triangles.
    """
    l = float(length)
    w = float(width)
    h = float(height)

    # 8 vertices
    # order: front face (s=0): 0..3, back face (s=l): 4..7
    verts = np.array([
        [-w/2,  h/2, 0.0],  # 0
        [ w/2,  h/2, 0.0],  # 1
        [ w/2, -h/2, 0.0],  # 2
        [-w/2, -h/2, 0.0],  # 3
        [-w/2,  h/2,  l ],  # 4
        [ w/2,  h/2,  l ],  # 5
        [ w/2, -h/2,  l ],  # 6
        [-w/2, -h/2,  l ],  # 7
    ], dtype=float)

    # 12 triangles (2 per face)
    # faces: front(0-1-2-3), back(4-5-6-7), left(0-3-7-4), right(1-2-6-5),
    #        top(0-1-5-4), bottom(3-2-6-7)
    tri = np.array([
        [0, 1, 2], [0, 2, 3],  # front
        [4, 6, 5], [4, 7, 6],  # back (note winding)
        [0, 3, 7], [0, 7, 4],  # left
        [1, 5, 6], [1, 6, 2],  # right
        [0, 4, 5], [0, 5, 1],  # top
        [3, 2, 6], [3, 6, 7],  # bottom
    ], dtype=int)

    i, j, k = tri[:, 0], tri[:, 1], tri[:, 2]
    return verts, i, j, k


def _cylinder_mesh_local(length: float, diameter: float, n: int = 28) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Cylinder aligned with local s-axis, s in [0, length].
    Returns:
      verts (N,3), i, j, k arrays for Mesh3d triangles.
    """
    l = float(length)
    r = 0.5 * float(diameter)
    n = int(max(8, n))

    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    ct, st = np.cos(theta), np.sin(theta)

    # rings
    front = np.stack([r * ct, r * st, np.zeros_like(ct)], axis=1)      # (n,3) s=0
    back  = np.stack([r * ct, r * st, np.full_like(ct, l)], axis=1)     # (n,3) s=l

    # centers for caps
    c_front = np.array([[0.0, 0.0, 0.0]], dtype=float)
    c_back  = np.array([[0.0, 0.0, l]], dtype=float)

    verts = np.vstack([front, back, c_front, c_back])  # (2n+2, 3)
    idx_front0 = 0
    idx_back0 = n
    idx_c_front = 2 * n
    idx_c_back = 2 * n + 1

    tris = []

    # side triangles: connect quads between rings
    for t in range(n):
        t2 = (t + 1) % n
        a = idx_front0 + t
        b = idx_front0 + t2
        c = idx_back0 + t2
        d = idx_back0 + t
        # two triangles (a-b-c) and (a-c-d)
        tris.append([a, b, c])
        tris.append([a, c, d])

    # front cap (winding so normal points outward-ish)
    for t in range(n):
        t2 = (t + 1) % n
        a = idx_front0 + t
        b = idx_front0 + t2
        tris.append([idx_c_front, b, a])

    # back cap
    for t in range(n):
        t2 = (t + 1) % n
        a = idx_back0 + t
        b = idx_back0 + t2
        tris.append([idx_c_back, a, b])

    tri = np.array(tris, dtype=int)
    i, j, k = tri[:, 0], tri[:, 1], tri[:, 2]
    return verts, i, j, k


def _transform_verts(r_start: np.ndarray, W: np.ndarray, verts_local: np.ndarray) -> np.ndarray:
    """
    Apply rigid transform: global = r_start + W @ local
    verts_local: (N,3)
    returns: (N,3)
    """
    r0 = np.asarray(r_start, dtype=float).reshape(3, 1)
    V = np.asarray(verts_local, dtype=float).T  # (3,N)
    G = r0 + W @ V
    return G.T


# -----------------------------
# Main plotting function
# -----------------------------
def plot_layout_3d(
    layout,
    show_orbit: bool = True,
    show_elements: bool = True,
    cylinder_points: int = 28,
    element_opacity: float = 0.95,
    orbit_width: int = 4,
    apply_tilt_to_body: bool = True,
    title: str = "Facility Layout (3D)",
):
    """
    Plot a MachineLayout in 3D using Plotly.

    apply_tilt_to_body:
      - If True: rotates element body meshes by element.tilt around local s.
        This is recommended for visualizing skew quads etc.
      - For bends, the survey already encodes the trajectory; the body tilt is still useful
        to show rotated yoke orientation.

    Returns:
      fig (plotly.graph_objects.Figure)
    """
    surveys = layout.survey()
    fig = go.Figure()

    # collect points for axis scaling
    all_x, all_y, all_z = [], [], []

    orbit_colors = ["black", "blue", "red", "green", "purple", "orange"]

    # -------- Orbits --------
    if show_orbit:
        for i, (line_name, data) in enumerate(surveys.items()):
            col = orbit_colors[i % len(orbit_colors)]

            X = [p["r_start"][0] for p in data]
            Y = [p["r_start"][1] for p in data]
            Z = [p["r_start"][2] for p in data]
            last = data[-1]
            X.append(last["r_end"][0])
            Y.append(last["r_end"][1])
            Z.append(last["r_end"][2])

            all_x.extend(X); all_y.extend(Y); all_z.extend(Z)

            fig.add_trace(go.Scatter3d(
                x=X, y=Y, z=Z,
                mode="lines",
                name=f"Orbit: {line_name}",
                line=dict(width=orbit_width, color=col),
                hoverinfo="name",
            ))

    # -------- Elements --------
    if show_elements:
        for line_name, data in surveys.items():
            for item in data:
                el = item.get("element", None)
                if el is None:
                    continue
                L = float(getattr(el, "l", 0.0))
                if L <= 0.0:
                    continue

                r0 = np.asarray(item["r_start"], dtype=float)
                r1 = np.asarray(item["r_end"], dtype=float)

                all_x.extend([r0[0], r1[0]])
                all_y.extend([r0[1], r1[1]])
                all_z.extend([r0[2], r1[2]])

                # Prefer W_start for geometry placement
                W_place = item.get("W_start", None)
                if W_place is None:
                    W_place = item["W"]
                W_place = np.asarray(W_place, dtype=float)

                # Visual params
                w = float(getattr(el, "width", 0.10))
                h = float(getattr(el, "height", 0.10))
                ang = float(getattr(el, "angle", 0.0))
                tilt = float(getattr(el, "tilt", 0.0))
                color = getattr(el, "color", None) or "lightgrey"

                cls_name = el.__class__.__name__
                is_cylinder = any(x in cls_name for x in ("Cavity", "Solenoid", "Drift"))

                # Build local mesh
                if is_cylinder:
                    verts_local, mi, mj, mk = _cylinder_mesh_local(L, diameter=w, n=cylinder_points)
                    # optionally rotate the body by tilt (purely visual)
                    if apply_tilt_to_body and tilt != 0.0:
                        T = _get_roll_matrix(tilt)
                        verts_local = (T @ verts_local.T).T
                else:
                    # For bends: approximate with a straight chord box and rotate it to chord direction
                    if ang != 0.0:
                        if abs(ang) < 1e-12:
                            l_chord = L
                        else:
                            l_chord = L * (np.sin(ang / 2.0) / (ang / 2.0))

                        verts_local, mi, mj, mk = _box_mesh_local(l_chord, w, h)

                        # Align box to chord direction: rot_ang = -angle/2 (your fix)
                        R_chord = _rot_y(-ang / 2.0)
                        verts_local = (R_chord @ verts_local.T).T

                        # additionally rotate body by tilt around s (visual yoke orientation)
                        if apply_tilt_to_body and tilt != 0.0:
                            T = _get_roll_matrix(tilt)
                            verts_local = (T @ verts_local.T).T
                    else:
                        verts_local, mi, mj, mk = _box_mesh_local(L, w, h)
                        if apply_tilt_to_body and tilt != 0.0:
                            T = _get_roll_matrix(tilt)
                            verts_local = (T @ verts_local.T).T

                # Transform to global
                verts_glob = _transform_verts(r0, W_place, verts_local)

                fig.add_trace(go.Mesh3d(
                    x=verts_glob[:, 0],
                    y=verts_glob[:, 1],
                    z=verts_glob[:, 2],
                    i=mi, j=mj, k=mk,
                    color=color,
                    opacity=element_opacity,
                    flatshading=True,
                    name=f"{line_name}.{getattr(el, 'id', cls_name)}",
                    hoverinfo="name",
                    lighting=dict(ambient=0.55, diffuse=0.75, specular=0.10, roughness=0.85, fresnel=0.02),
                ))

    # -------- Make X and Y scales equal (fix “round tubes look elliptical”) --------
    if not all_x:
        all_x = [0.0]; all_y = [0.0]; all_z = [0.0]

    xmin, xmax = float(min(all_x)), float(max(all_x))
    ymin, ymax = float(min(all_y)), float(max(all_y))
    zmin, zmax = float(min(all_z)), float(max(all_z))

    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)

    rx = max(xmax - xmin, 1e-6)
    ry = max(ymax - ymin, 1e-6)
    rz = max(zmax - zmin, 1e-6)

    # enforce minimum transverse span for “almost perfectly straight” machines
    span_xy = max(rx, ry, 0.5)
    pad_xy = 0.05 * span_xy
    half_xy = 0.5 * span_xy + pad_xy

    x_range = [cx - half_xy, cx + half_xy]
    y_range = [cy - half_xy, cy + half_xy]

    pad_z = 0.02 * rz
    z_range = [zmin - pad_z, zmax + pad_z]
    rz_padded = (z_range[1] - z_range[0])

    aspect_z = max(rz_padded / span_xy, 0.2)

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis=dict(title="X [m]", range=x_range),
            yaxis=dict(title="Y [m]", range=y_range),
            zaxis=dict(title="Z [m]", range=z_range),
            aspectmode="manual",
            aspectratio=dict(x=1.0, y=1.0, z=aspect_z),
            camera=dict(
                up=dict(x=0, y=1, z=0),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=1.6, y=0.9, z=0.55),
            ),
        ),
        margin=dict(l=0, r=0, b=0, t=45),
        legend=dict(itemsizing="constant"),
    )

    return fig


# -----------------------------
# Simple demo usage (optional)
# -----------------------------
if __name__ == "__main__":
    # Minimal self-test if you run this module directly (requires Ocelot objects)
    from ocelot import Drift, Quadrupole, SBend
    from ocelot.cpbd.magnetic_lattice import MagneticLattice

    # import your MachineLayout from wherever you placed it:
    # from ocelot.cpbd.machine_layout import MachineLayout
    from ocelot.cpbd.layout import MachineLayout  # adjust if needed

    main = MagneticLattice([
        Drift(l=5, eid="d1"),
        Quadrupole(l=0.5, eid="q1", tilt=np.pi/4, width=0.20, height=0.20, color="orange"),
        Drift(l=1, eid="d2"),
        SBend(l=2, angle=0.25, eid="b1", width=0.30, height=0.30, color="steelblue"),
        Drift(l=3, eid="d3"),
        Quadrupole(l=0.5, eid="q2", width=0.20, height=0.20, color="orange"),
    ])

    dump = MagneticLattice([
        Drift(l=1, eid="dd1"),
        SBend(l=2, angle=0.45, eid="db1", tilt=np.pi/2, width=0.25, height=0.25, color="crimson"),
        Drift(l=4, eid="dd2"),
    ])

    facility = MachineLayout()
    facility.add_line("Main", main)
    facility.add_line("Dump", dump, parent_name="Main", anchor_element_id="b1")

    fig = plot_layout_3d(facility, show_orbit=True, show_elements=True)
    fig.show()