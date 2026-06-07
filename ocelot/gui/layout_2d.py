import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


def _rect_corners_2d(p0, p1, thickness):
    """
    Build 4 corners of a 2D rectangle centered on segment p0->p1.
    p0, p1: (2,) arrays
    thickness: full thickness (not radius)
    Returns (4,2) array of corners.
    """
    v = p1 - p0
    L = float(np.hypot(v[0], v[1]))
    if L < 1e-12:
        # Degenerate element: return a tiny square
        eps = 0.5 * max(thickness, 1e-3)
        return np.array([
            [p0[0] - eps, p0[1] - eps],
            [p0[0] + eps, p0[1] - eps],
            [p0[0] + eps, p0[1] + eps],
            [p0[0] - eps, p0[1] + eps],
        ])

    t = v / L
    n = np.array([-t[1], t[0]])  # left normal
    half = 0.5 * thickness

    c0 = p0 + n * half
    c1 = p1 + n * half
    c2 = p1 - n * half
    c3 = p0 - n * half
    return np.vstack([c0, c1, c2, c3])


def plot_layout_2d(
    layout,
    show_orbit=True,
    show_elements=True,
    element_alpha=0.35,
    orbit_lw=2.0,
    element_edge_lw=0.6,
    equal_aspect=True,
    title="Facility Layout (2D: Top + Side)"
):
    """
    Plots the facility in 2D using Matplotlib:
      - Top view:  X vs Z
      - Side view: Y vs Z

    Requires layout.survey() to have been implemented.
    Uses:
      - item['r_start'], item['r_end'] (3-vectors)
      - item['element'] with .l, .width, .height, .color
      - item.get('W_start') is optional (not strictly required for 2D projection)

    Returns: (fig, (ax_top, ax_side))
    """
    surveys = layout.survey()  # ensure fresh

    fig, (ax_top, ax_side) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle(title)

    # Color cycle for orbits (independent of element colors)
    orbit_colors = ["black", "tab:blue", "tab:red", "tab:green", "tab:purple", "tab:orange"]

    # --- A) Orbits ---
    if show_orbit:
        for i, (line_name, data) in enumerate(surveys.items()):
            col = orbit_colors[i % len(orbit_colors)]

            # polyline of starts + last end
            X = [p["r_start"][0] for p in data]
            Y = [p["r_start"][1] for p in data]
            Z = [p["r_start"][2] for p in data]

            last = data[-1]
            X.append(last["r_end"][0])
            Y.append(last["r_end"][1])
            Z.append(last["r_end"][2])

            ax_top.plot(Z, X, lw=orbit_lw, color=col, label=f"Orbit: {line_name}")
            ax_side.plot(Z, Y, lw=orbit_lw, color=col, label=f"Orbit: {line_name}")

    # --- B) Elements as rectangles projected to each view ---
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

                # Thickness in each view
                w = float(getattr(el, "width", 0.05))   # used in top view as thickness in X
                h = float(getattr(el, "height", 0.05))  # used in side view as thickness in Y

                face = getattr(el, "color", None) or "lightgray"

                # TOP view rectangle in (Z, X)
                p0_top = np.array([r0[2], r0[0]])  # (Z, X)
                p1_top = np.array([r1[2], r1[0]])  # (Z, X)
                corners_top = _rect_corners_2d(p0_top, p1_top, thickness=w)

                ax_top.add_patch(
                    Polygon(
                        corners_top,
                        closed=True,
                        facecolor=face,
                        edgecolor="k",
                        linewidth=element_edge_lw,
                        alpha=element_alpha,
                    )
                )

                # SIDE view rectangle in (Z, Y)
                p0_side = np.array([r0[2], r0[1]])  # (Z, Y)
                p1_side = np.array([r1[2], r1[1]])  # (Z, Y)
                corners_side = _rect_corners_2d(p0_side, p1_side, thickness=h)

                ax_side.add_patch(
                    Polygon(
                        corners_side,
                        closed=True,
                        facecolor=face,
                        edgecolor="k",
                        linewidth=element_edge_lw,
                        alpha=element_alpha,
                    )
                )

                # Optional: label (can get crowded; comment out if not needed)
                # mid = 0.5 * (r0 + r1)
                # ax_top.text(mid[2], mid[0], el.id, fontsize=7, rotation=0)

    # --- Axes cosmetics ---
    ax_top.set_ylabel("X [m] (Top view)")
    ax_side.set_ylabel("Y [m] (Side view)")
    ax_side.set_xlabel("Z [m]")

    ax_top.grid(True, alpha=0.25)
    ax_side.grid(True, alpha=0.25)

    # Legends (avoid duplicates)
    if show_orbit:
        ax_top.legend(loc="best", fontsize=9)

    if equal_aspect:
        # For each view, keep equal data scale in both directions (optional)
        ax_top.set_aspect("equal", adjustable="datalim")
        ax_side.set_aspect("equal", adjustable="datalim")

    plt.tight_layout()
    return fig, (ax_top, ax_side)

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

    plot_layout_2d(facility, show_orbit=True, show_elements=True)
    plt.show()