__author__ = 'Sergey Tomin'

import copy
import sys
from argparse import ArgumentParser
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np

from ocelot import *
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaGlobalTM, RungeKuttaOcelotTM


COORD_LABELS = ("x", "px", "y", "py", "tau", "p")


def make_probe_coordinates():
    """
    Small offsets in Ocelot coordinates.

    Particle 0 is the zero particle.  It is useful for looking at what happens
    to the reference orbit in global RK tracking.
    """
    return np.array([
        [0.0,      0.0,      0.0,      0.0,      0.0, 0.0],
        [1.0e-4,  0.0,      0.0,      0.0,      0.0, 0.0],
        [-1.0e-4, 0.0,      0.0,      0.0,      0.0, 0.0],
        [0.0,      2.0e-5,  0.0,      0.0,      0.0, 0.0],
        [0.0,     -2.0e-5,  0.0,      0.0,      0.0, 0.0],
        [0.0,      0.0,      1.0e-4,  0.0,      0.0, 0.0],
        [0.0,      0.0,     -1.0e-4,  0.0,      0.0, 0.0],
        [0.0,      0.0,      0.0,      2.0e-5,  0.0, 0.0],
        [0.0,      0.0,      0.0,     -2.0e-5,  0.0, 0.0],
    ])


def make_parray(coordinates, energy):
    parray = ParticleArray(n=len(coordinates))
    parray.E = energy
    parray.rparticles[:] = coordinates.T
    return parray


def track_element(element, tm_class, parray):
    tracked = copy.deepcopy(parray)
    element.set_tm(tm_class)
    for tm in element.tms:
        tm.apply(tracked)
    return tracked.rparticles.copy()


def compare_element(name, element_factory, coordinates, energy):
    parray = make_parray(coordinates, energy)

    second = track_element(element_factory(), SecondTM, parray)
    global_rk = track_element(element_factory(), RungeKuttaGlobalTM, parray)
    legacy_rk = track_element(element_factory(), RungeKuttaTM, parray)
    ocelot_rk = track_element(element_factory(), RungeKuttaOcelotTM, parray)

    second_probe = second[:, 1:]
    global_probe = global_rk[:, 1:]
    ocelot_probe = ocelot_rk[:, 1:]

    return {
        "name": name,
        "global_ref": global_rk[:, 0],
        "ocelot_ref": ocelot_rk[:, 0],
        "legacy_matches_global": np.allclose(legacy_rk, global_rk, rtol=0.0, atol=0.0),
        "global_diff": global_probe - second_probe,
        "global_minus_ref_diff": global_probe - global_rk[:, [0]] - second_probe,
        "ocelot_diff": ocelot_probe - second_probe,
    }


def print_summary(results):
    print("SecondTM tracks in Ocelot coordinates: offsets relative to the design orbit.")
    print("RungeKuttaGlobalTM tracks in the fixed field frame: useful for trajectories, SR, CSR.")
    print("RungeKuttaTM is kept as the legacy name for RungeKuttaGlobalTM.")
    print("RungeKuttaOcelotTM tracks through the same field, then converts back to Ocelot coordinates.\n")

    for result in results:
        global_max = np.max(np.abs(result["global_diff"]), axis=1)
        global_minus_ref_max = np.max(np.abs(result["global_minus_ref_diff"]), axis=1)
        ocelot_max = np.max(np.abs(result["ocelot_diff"]), axis=1)

        print(result["name"])
        print("  Global RK reference particle at exit:")
        for label, value in zip(COORD_LABELS, result["global_ref"]):
            print(f"    {label:>3s}: {value: .8e}")
        print("  Ocelot RK reference particle at exit:")
        for label, value in zip(COORD_LABELS, result["ocelot_ref"]):
            print(f"    {label:>3s}: {value: .8e}")
        print(f"  Legacy RungeKuttaTM matches RungeKuttaGlobalTM: {result['legacy_matches_global']}")
        print("  max |Global RK - SecondTM|:")
        print("   ", np.array2string(global_max, precision=3, suppress_small=False))
        print("  max |(Global RK - global reference) - SecondTM|:")
        print("   ", np.array2string(global_minus_ref_max, precision=3, suppress_small=False))
        print("  max |Ocelot RK - SecondTM|:")
        print("   ", np.array2string(ocelot_max, precision=3, suppress_small=False))
        print()


def plot_summary(results):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    x = np.arange(len(COORD_LABELS))
    width = 0.28

    fig, axes = plt.subplots(len(results), 1, figsize=(9, 8), sharex=True)
    if len(results) == 1:
        axes = [axes]

    for ax, result in zip(axes, results):
        global_max = np.max(np.abs(result["global_diff"]), axis=1)
        global_minus_ref_max = np.max(np.abs(result["global_minus_ref_diff"]), axis=1)
        ocelot_max = np.max(np.abs(result["ocelot_diff"]), axis=1)

        ax.bar(x - width, global_max, width=width, label="Global RK")
        ax.bar(x, global_minus_ref_max, width=width, label="Global RK - reference")
        ax.bar(x + width, ocelot_max, width=width, label="Ocelot RK")
        ax.set_yscale("log")
        ax.set_ylabel(result["name"])
        ax.grid(True, which="both", axis="y", alpha=0.3)
        ax.legend(loc="best")

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(COORD_LABELS)
    axes[-1].set_xlabel("coordinate")
    fig.suptitle("Runge-Kutta tracking compared with SecondTM")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = ArgumentParser(description="Compare matrix tracking with global-frame and Ocelot-coordinate RK tracking.")
    parser.add_argument("--plot", action="store_true", help="Show a matplotlib summary plot.")
    args = parser.parse_args()

    energy = 1.0
    coordinates = make_probe_coordinates()

    elements = (
        ("Bend", lambda: Bend(l=0.5, angle=0.01, e1=0.0, e2=0.0, eid="B")),
        ("Quadrupole", lambda: Quadrupole(l=0.5, k1=0.7, eid="Q")),
        ("Sextupole", lambda: Sextupole(l=0.2, k2=10.0, eid="SX")),
    )

    results = [
        compare_element(name, element_factory, coordinates, energy)
        for name, element_factory in elements
    ]

    print_summary(results)
    if args.plot:
        plot_summary(results)
