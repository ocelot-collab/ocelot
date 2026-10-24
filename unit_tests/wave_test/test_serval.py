"""Small-grid checks for the SERVAL generator's public calling convention."""

import numpy as np
import pytest

from ocelot.common.globals import h_eV_s, speed_of_light
from ocelot.optics import wave


@pytest.mark.parametrize("plot_options", [{}, {"showfig": False}, {"showfig": True}])
def test_serval_grid_and_optional_plot(monkeypatch, plot_options):
    from ocelot.gui import dfl_plot

    plotted = []
    monkeypatch.setattr(dfl_plot, "plot_dfl", lambda field, **kwargs: plotted.append(field))
    monkeypatch.setattr(wave, "nthread", 1)
    field = wave.dfl_gen_undulator_serval(
        E_ph=2167, L_w=3.6, shape=(15, 17, 2), dgrid=(400e-6, 300e-6, 1),
        sig_x=38e-6, sig_y=4.68e-6, sig_xp=25e-6, sig_yp=20e-6,
        seed=1, **plot_options,
    )
    assert field.fld.shape == (2, 17, 15)
    assert field.dx == pytest.approx(400e-6 / 15)
    assert field.dy == pytest.approx(300e-6 / 17)
    assert field.dz == pytest.approx(0.5)
    assert field.xlamds == pytest.approx(h_eV_s * speed_of_light / 2167)
    assert np.isfinite(field.fld).all()
    assert np.any(field.fld != 0)
    assert len(plotted) == int(plot_options.get("showfig", False))
