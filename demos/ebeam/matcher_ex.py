"""Ultimate matcher example.

This demo shows three key matcher capabilities in one script:

1) Custom variable (`Vary`): one chicane knob controls multiple lattice values
   consistently (4 dipole angles, edge angles, and shoulder drifts).
2) Custom objective (`Objective`): optimize peak current using full particle
   tracking inside matcher residual evaluation.
3) Standard target (`target_twiss`): constrain final beam energy.

Run:
    python demos/ebeam/matcher_ex.py

Output:
- Console summary with optimized knobs and achieved peak current.
- Figure file `matcher_ex_before_after.png` with before/after current profile
  and longitudinal phase space.
"""

from __future__ import annotations

import copy
import os
import sys
from pathlib import Path

# Ensure local checkout is imported, not an older site-packages installation.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np

# Use a non-interactive backend only for truly headless Linux sessions.
if (
    "MPLBACKEND" not in os.environ
    and sys.platform.startswith("linux")
    and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
):
    import matplotlib

    matplotlib.use("Agg")

import matplotlib.pyplot as plt

from ocelot import (
    Cavity,
    Drift,
    MagneticLattice,
    Marker,
    SBend,
    Twiss,
    generate_parray,
    track,
)
from ocelot.cpbd.matcher import MatchProblem, Objective, Vary


# -----------------------------------------------------------------------------
# User configuration
# -----------------------------------------------------------------------------
TARGET_PEAK_CURRENT_A = 300.0
TARGET_END_ENERGY_GEV = 0.13

THETA0 = 0.10
THETA_LIMITS = (0.01, 0.30)

NPARTICLES = 10000
TOTAL_CHARGE_C = 250e-12
SIGMA_TAU_M = 1e-3
SIGMA_P = 5e-4

RNG_SEED = 7
MAX_ITER = 300



class PeakCurrentObjective(Objective):
    """Tracking-based objective for peak current matching."""

    def __init__(self, parray_template, target_current: float, num_bins: int = 200, **kwargs):
        super().__init__(**kwargs)
        self.parray_template = parray_template
        self.target_current = float(target_current)
        self.num_bins = int(num_bins)

    def residuals(self, state):
        # Track a fresh copy each evaluation so residuals always use the same
        # initial particle distribution.
        parray = copy.deepcopy(self.parray_template)
        _tws_track, tracked = track(state.lat, parray, print_progress=False)
        i_peak = peak_current(tracked, num_bins=self.num_bins)
        return np.array([(i_peak - self.target_current) / self.target_current], dtype=float)


# -----------------------------------------------------------------------------
# Main demo
# -----------------------------------------------------------------------------
def main() -> None:
    np.random.seed(RNG_SEED)

    # Build a C-shape chicane section.
    z_distance = 1.0
    d_ref = Drift(l=0.5)
    d12 = Drift(l=shoulders_from_theta(THETA0, z_distance))
    d34 = Drift(l=shoulders_from_theta(THETA0, z_distance))

    b1 = SBend(l=0.2, angle=-THETA0, e1=0.0, e2=-THETA0, tilt=0.0, fint=0.0, eid="B1")
    b2 = SBend(l=0.2, angle=+THETA0, e1=THETA0, e2=0.0, tilt=0.0, fint=0.0, eid="B2")
    b3 = SBend(l=0.2, angle=+THETA0, e1=0.0, e2=THETA0, tilt=0.0, fint=0.0, eid="B3")
    b4 = SBend(l=0.2, angle=-THETA0, e1=-THETA0, e2=0.0, tilt=0.0, fint=0.0, eid="B4")
    chicane = (b1, d12, b2, d_ref, b3, d34, b4)

    # RF section + end marker.
    c11 = Cavity(l=0.5, v=0.15, phi=0.0, freq=1.3e9, eid="C11")
    c13 = Cavity(l=0.5, v=0.0166, phi=180.0, freq=3.9e9, eid="C13")
    end = Marker(eid="END")

    cell = [d_ref, c11, d_ref, c13, d_ref, *chicane, d_ref, end]
    lat = MagneticLattice(cell)

    tw0 = Twiss(beta_x=10.0, beta_y=10.0, emit_xn=1e-6, emit_yn=1e-6, E=0.005)
    parray_init = generate_parray(
        tws=tw0,
        nparticles=NPARTICLES,
        charge=TOTAL_CHARGE_C,
        chirp=0.0,
        sigma_p=SIGMA_P,
        sigma_tau=SIGMA_TAU_M,
    )

    # Track once before matching for reference.
    _tb, tracked_before = track(lat, copy.deepcopy(parray_init), print_progress=False)
    i_before = peak_current(tracked_before)

    # Build matcher problem.
    problem = MatchProblem(lat, tw0)

    # Custom composite variable: one theta knob drives whole chicane geometry.
    def get_theta() -> float:
        return float(b2.angle)

    def set_theta(theta: float) -> None:
        # 1) Linked dipole angles.
        b1.angle = -theta
        b2.angle = +theta
        b3.angle = +theta
        b4.angle = -theta

        # 2) Linked edge angles.
        b1.e1, b1.e2 = 0.0, -theta
        b2.e1, b2.e2 = +theta, 0.0
        b3.e1, b3.e2 = 0.0, +theta
        b4.e1, b4.e2 = -theta, 0.0

        # 3) Recompute shoulder drifts from geometry.
        l_shoulder = shoulders_from_theta(theta, z_distance)
        d12.l = l_shoulder
        d34.l = l_shoulder

    problem.add_variable(
        Vary(
            name="CHICANE_theta",
            getter=get_theta,
            setter=set_theta,
            limits=THETA_LIMITS,
        )
    )

    # Standard element variables.
    problem.vary_element(c11, "v", limits=(0.10, 0.20), name="C11_v")
    problem.vary_element(c11, "phi", limits=(-40.0, 40.0), name="C11_phi")
    problem.vary_element(c13, "v", limits=(0.01, 0.025), name="C13_v")
    problem.vary_element(c13, "phi", limits=(90.0, 270.0), name="C13_phi")

    # Custom tracking objective: target peak current.
    problem.add_objective(
        PeakCurrentObjective(
            parray_template=parray_init,
            target_current=TARGET_PEAK_CURRENT_A,
            name="peak_current_obj",
            weight=1e1,
        )
    )

    # Standard physics constraint: final energy at end marker.
    problem.target_twiss(end, "E", value=TARGET_END_ENERGY_GEV, weight=1e6)

    # Solve.
    result = problem.solve(solver="ls_trf", max_iter=MAX_ITER)

    # Track after matching.
    _ta, tracked_after = track(lat, copy.deepcopy(parray_init), print_progress=False)
    i_after = peak_current(tracked_after)

    # Final energy (from matched Twiss state used by matcher).
    _merit, _trep, _orep, state = problem.evaluate()
    end_energy = float(state.twiss_end.E)

    print("\n=== Matcher Ultimate Example Summary ===")
    print(f"Solve success: {result.success}")
    print(f"Target peak current [A]: {TARGET_PEAK_CURRENT_A:.3f}")
    print(f"Initial peak current [A]: {i_before:.3f}")
    print(f"Final peak current [A]: {i_after:.3f}")
    print(f"Final end energy [GeV]: {end_energy:.6f} (target {TARGET_END_ENERGY_GEV:.6f})")
    print("Optimized variables:")
    for k, v in result.variables.items():
        print(f"  {k:>14s} = {v:.8g}")

    out_png = Path(__file__).with_name("matcher_ex_before_after.png")
    make_before_after_plot(tracked_before, tracked_after, out_png)
    print(f"Saved plot: {out_png}")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def shoulders_from_theta(theta: float, z_distance: float) -> float:
    return z_distance / np.cos(theta)


def peak_current(parray, num_bins: int = 200) -> float:
    return float(np.max(parray.I(num_bins=num_bins)[:, 1]))


def make_before_after_plot(parray_before, parray_after, output_path: Path) -> None:
    cur_before = parray_before.I(num_bins=200)
    cur_after = parray_after.I(num_bins=200)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    axes[0, 0].plot(cur_before[:, 0] * 1e3, cur_before[:, 1], lw=1.5)
    axes[0, 0].set_title("Current Profile (Before)")
    axes[0, 0].set_xlabel("tau [mm]")
    axes[0, 0].set_ylabel("I [A]")

    axes[0, 1].plot(cur_after[:, 0] * 1e3, cur_after[:, 1], lw=1.5, color="tab:orange")
    axes[0, 1].set_title("Current Profile (After)")
    axes[0, 1].set_xlabel("tau [mm]")
    axes[0, 1].set_ylabel("I [A]")

    axes[1, 0].scatter(parray_before.tau() * 1e3, parray_before.p(), s=2, alpha=0.25)
    axes[1, 0].set_title("Longitudinal Phase Space (Before)")
    axes[1, 0].set_xlabel("tau [mm]")
    axes[1, 0].set_ylabel("dp/p")

    axes[1, 1].scatter(parray_after.tau() * 1e3, parray_after.p(), s=2, alpha=0.25, color="tab:orange")
    axes[1, 1].set_title("Longitudinal Phase Space (After)")
    axes[1, 1].set_xlabel("tau [mm]")
    axes[1, 1].set_ylabel("dp/p")

    fig.suptitle("Matcher Ultimate Example: Before vs After", fontsize=13)
    fig.savefig(output_path, dpi=150)

    # `Agg` (and similar) cannot open a GUI window; keep interactive show only
    # when a GUI backend is active.
    if "agg" in plt.get_backend().lower():
        plt.close(fig)
    else:
        plt.show()
        plt.close(fig)


if __name__ == "__main__":
    main()
