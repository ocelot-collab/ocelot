"""Survey and 2D plotting must remain available without the optional 3D backend."""

import os
from pathlib import Path
import subprocess
import sys
import textwrap


def test_layout_without_plotly(tmp_path):
    code = textwrap.dedent('''
        import sys
        from importlib.abc import MetaPathFinder

        class WithoutPlotly(MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname == "plotly" or fullname.startswith("plotly."):
                    raise ModuleNotFoundError("Plotly is not installed", name=fullname)

        sys.meta_path.insert(0, WithoutPlotly())

        import matplotlib.pyplot as plt
        from ocelot.cpbd.elements import Drift
        from ocelot.cpbd.layout import MachineLayout
        from ocelot.cpbd.magnetic_lattice import MagneticLattice
        from ocelot.gui.layout_2d import plot_layout_2d
        from ocelot.gui.layout_3d import plot_layout_3d

        layout = MachineLayout()
        layout.add_line("line", MagneticLattice([Drift(l=2)]))
        assert layout.survey()["line"][-1]["Z"] == 2
        fig, axes = plot_layout_2d(layout)
        fig.canvas.draw()
        assert len(axes) == 2
        plt.close(fig)

        try:
            plot_layout_3d(layout)
        except ImportError as exc:
            assert 'pip install "ocelot-collab[plotly]"' in str(exc)
        else:
            raise AssertionError("3D plotting should require Plotly")
    ''')
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=Path(__file__).resolve().parents[2],
        env={**os.environ, "MPLBACKEND": "Agg", "MPLCONFIGDIR": str(tmp_path)},
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_layout_with_plotly():
    import pytest

    pytest.importorskip("plotly")
    from ocelot.cpbd.elements import Drift, Quadrupole, SBend
    from ocelot.cpbd.layout import MachineLayout
    from ocelot.cpbd.magnetic_lattice import MagneticLattice
    from ocelot.gui.layout_3d import plot_layout_3d

    layout = MachineLayout()
    layout.add_line("line", MagneticLattice([
        Drift(l=1), Quadrupole(l=0.2, k1=1), SBend(l=0.5, angle=0.1),
    ]))
    fig = plot_layout_3d(layout)
    assert any(trace.type == "scatter3d" for trace in fig.data)
    assert any(trace.type == "mesh3d" for trace in fig.data)
    assert '"data"' in fig.to_json()
