"""Measure warm CSR kicks at fixed resolution, optionally against another csr.py.

From the repository root::

    python -m benchmarks.benchmark_csr --output /tmp/csr-timing.json
    python -m benchmarks.benchmark_csr --reference /tmp/csr-before.py
    python -m benchmarks.benchmark_csr --profile-dir /tmp/csr-profiles

The timed operation is one 5 cm kick inside a bend, not full lattice tracking.
Compilation, trajectory preparation, particle generation and momentum resets are
outside the timed samples. Reference/current calls alternate to reduce timing
bias. All cases use identical physical and numerical settings.
"""

import argparse
import cProfile
import importlib.util
import json
import platform
import sys
from pathlib import Path
from statistics import median
from time import perf_counter

import numpy as np

from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.csr import CSR, nb_flag, ne_flag
from ocelot.cpbd.elements import Bend, Drift
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.navi import Navigator


def make_process(cls, lattice, args):
    process = cls(n_bin=args.bins, m_bin=5, sigma_min=args.sigma_min,
                  ip_method=2, sp=.5, traj_step=.0002, apply_step=.0005)
    # Only reference modules from before the rename use the old name.
    if hasattr(process, 'mesh_step_quantum'):
        process.mesh_step_quantum = 0
    else:
        process.step_unit = 0
    Navigator(lattice).add_physics_proc(process, lattice.sequence[0], lattice.sequence[-1])
    process.z0 = .4
    return process


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--sizes', type=int, nargs='+', default=[10000, 100000, 1000000])
    parser.add_argument('--charge-binning-fraction',
                        type=float, nargs='+', default=[0., .5],
                        help='0: length bins, 1: charge bins; intermediate values blend both')
    parser.add_argument('--bins', type=int, default=100)
    parser.add_argument('--sigma-min', type=float, default=2e-5)
    parser.add_argument('--repeats', type=int, default=9)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--profile-dir', type=Path, help='save a warm current-kick cProfile per case')
    args = parser.parse_args()
    if args.repeats < 1 or min(args.sizes) < 2:
        parser.error('repeats must be positive and particle counts at least two')
    if args.profile_dir:
        args.profile_dir.mkdir(parents=True, exist_ok=True)

    lattice = MagneticLattice([Drift(l=.1, eid='bench_in'),
                               Bend(l=.5, angle=.1, eid='bench_bend'),
                               Drift(l=1., eid='bench_out')])
    classes = {'current': CSR}
    if args.reference:
        spec = importlib.util.spec_from_file_location('csr_benchmark_reference', args.reference)
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        classes['reference'] = module.CSR
    processes = {}
    setup = {}
    for name, cls in classes.items():
        start = perf_counter()
        processes[name] = make_process(cls, lattice, args)
        setup[name] = perf_counter() - start

    report = {'platform': platform.platform(), 'python': platform.python_version(),
              'numpy': np.__version__, 'numba_enabled': nb_flag, 'numexpr_enabled': ne_flag,
              'bins': args.bins, 'm_bin': 5, 'sigma_min_m': args.sigma_min,
              'ip_method': 2, 'sp': .5, 'mesh_step_quantum_m': 0,
              'traj_step_m': .0002, 'apply_step_m': .0005, 'kick_length_m': .05,
              'energy_GeV': .13, 'sigma_z_m': 2e-4, 'charge_C': 1e-9,
              'setup_seconds': setup, 'cases': []}
    for count in args.sizes:
        rng = np.random.default_rng(711)
        particles = ParticleArray(count)
        particles.E = .13
        particles.rparticles[4] = rng.normal(0, 2e-4, count)
        particles.q_array[:] = rng.uniform(.1, 3, count)
        particles.q_array *= 1e-9 / particles.q_array.sum()
        for x_qbin in args.charge_binning_fraction:
            times = {name: [] for name in processes}
            kicks = {}
            first_call = {}
            for name, process in processes.items():
                if hasattr(process, 'charge_binning_fraction'):
                    process.charge_binning_fraction = x_qbin
                else:
                    process.x_qbin = x_qbin
                particles.p().fill(0)
                start = perf_counter()
                process.apply(particles, .05)
                first_call[name] = perf_counter() - start
                kicks[name] = particles.p().copy()
            for repeat in range(args.repeats):
                names = list(processes)
                if repeat % 2:
                    names.reverse()
                for name in names:
                    particles.p().fill(0)
                    start = perf_counter()
                    processes[name].apply(particles, .05)
                    times[name].append(perf_counter() - start)
            row = {'particles': count, 'charge_binning_fraction': x_qbin,
                   'first_call_seconds': first_call,
                   'median_seconds': {name: median(values) for name, values in times.items()},
                   'samples_seconds': times}
            if 'reference' in kicks:
                difference = kicks['current'] - kicks['reference']
                row['relative_kick_l2_difference'] = float(
                    np.linalg.norm(difference) / np.linalg.norm(kicks['reference']))
                row['max_abs_dp_difference'] = float(np.max(np.abs(difference)))
                row['speedup'] = row['median_seconds']['reference'] / row['median_seconds']['current']
                np.testing.assert_allclose(kicks['current'], kicks['reference'], rtol=1e-11, atol=1e-16)
            print(json.dumps(row), flush=True)
            report['cases'].append(row)
            if args.profile_dir:
                particles.p().fill(0)
                profile = cProfile.Profile()
                profile.runcall(processes['current'].apply, particles, .05)
                profile.dump_stats(args.profile_dir / f'csr-{count}-x{x_qbin}.prof')
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
