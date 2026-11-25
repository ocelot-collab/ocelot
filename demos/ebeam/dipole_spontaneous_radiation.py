"""Track a monoenergetic bunch through a dipole with the undulator parameter K set to zero."""

import matplotlib.pyplot as plt
import numpy as np

import ocelot as ocl
from ocelot.common.globals import lambda_C_r, m_e_GeV, ro_e
from ocelot.gui.accelerator import plot_API


def main():
    entrance = ocl.Drift(l=0.5, eid="Entrance drift")
    bend = ocl.SBend(l=1.0, angle=0.01, eid="Dipole")
    exit_drift = ocl.Drift(l=0.5, eid="Exit drift")
    lattice = ocl.MagneticLattice((entrance, bend, exit_drift))

    # The process runs from the bend entrance to the entrance of the exit drift.
    # Only radius sets the dipole strength. K=0 is also the default.
    radius = bend.l / bend.angle
    radiation = ocl.SpontanRadEffects(type="dipole", radius=radius, K=0)
    navi = ocl.Navigator(lattice, unit_step=0.05)
    navi.add_physics_proc(radiation, bend, exit_drift)

    np.random.seed(17)
    beam = ocl.generate_parray(energy=10., sigma_p=0., chirp=0., nparticles=40_000)
    tws, _ = ocl.track(lattice, beam, navi, print_progress=False)
    s = np.array([t.s for t in tws])
    # The reference energy stays fixed. Convert the tracked energy moments to keV.
    mean_loss = -beam.p0c * 1e6 * np.array([t.p for t in tws])
    rms_spread = beam.p0c * 1e6 * np.sqrt([t.pp for t in tws])

    # Uniform-dipole predictions at the entrance energy (small fractional loss).
    s_curve = np.linspace(0., lattice.totalLen, 401)
    magnetic_length = np.clip(s_curve - entrance.l, 0., bend.l)
    gamma = beam.E / m_e_GeV
    loss_per_m = (2 / 3) * ro_e * gamma**3 / radius**2
    variance_per_m = 55 / (24 * np.sqrt(3)) * ro_e * lambda_C_r * gamma**5 / abs(radius)**3
    predicted_loss = beam.E * 1e6 * loss_per_m * magnetic_length
    predicted_spread = beam.E * 1e6 * np.sqrt(variance_per_m * magnetic_length)

    fig, (ax_loss, ax_spread) = plot_API(lattice, add_extra_subplot=True, legend=False,
                                         figsize=(8, 6), fig_name="Dipole radiation")
    fig.suptitle("Spontaneous radiation in a 10 GeV beam  ·  1 m dipole, K = 0")
    for ax, tracked, predicted, label in (
        (ax_loss, mean_loss, predicted_loss, "Mean energy loss [keV]"),
        (ax_spread, rms_spread, predicted_spread, "Rms energy spread [keV]"),
    ):
        ax.plot(s_curve, predicted, color="C0", label="Uniform-dipole theory")
        ax.plot(s, tracked, "o", color="C1", markersize=4, label="Tracked bunch")
        ax.set_ylabel(label)
    ax_loss.legend(loc="upper left")
    ax_spread.set_xlim(0., lattice.totalLen)
    print(f"Exit mean loss: {mean_loss[-1]:.2f} keV (theory {predicted_loss[-1]:.2f} keV)")
    print(f"Exit rms spread: {rms_spread[-1]:.2f} keV (theory {predicted_spread[-1]:.2f} keV)")
    plt.show()


if __name__ == "__main__":
    main()
