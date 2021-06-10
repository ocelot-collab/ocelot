from pathlib import Path
import pytest

import numpy as np
from numpy.linalg import norm
import easygdf

from ocelot.adaptors import gpt as gpt_adaptor
from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.beam import ParticleArray


# Toy data for instantiating a Screen
TEST_GDF_FILE_PATH = Path(__file__).parent / "gpt-2-particles-1-screen-1-tout.gdf"


@pytest.fixture
def example_gdf():
    """tmp path to written pmd h5 file defined from above PARTICLE_GROUP_DATA"""
    gdf = easygdf.load_screens_touts(str(TEST_GDF_FILE_PATH))
    yield gdf


@pytest.fixture
def example_screen(example_gdf):
    yield example_gdf["screens"][0]


@pytest.fixture
def example_tout(example_gdf):
    yield example_gdf["touts"][0]


def compare_gdf_dict_with_parray(gdfmap, parray):
    energy = gdfmap["G"] * m_e_GeV
    momentum = (energy ** 2 - m_e_GeV ** 2) ** 0.5
    reference_energy = np.mean(energy)
    p0 = np.mean(momentum)

    if "time" in gdfmap: # tout
        tau = np.mean(gdfmap["z"]) - gdfmap["z"]
    else: # screen
        dt = gdfmap["t"] - np.mean(gdfmap["t"])
        tau = speed_of_light * dt

    assert parray.E == reference_energy
    np.testing.assert_allclose(parray.x(), gdfmap["x"] - np.mean(gdfmap["x"]))
    np.testing.assert_allclose(parray.px(), gdfmap["Bx"] * momentum / p0)
    np.testing.assert_allclose(parray.y(), gdfmap["y"] - np.mean(gdfmap["y"]))
    np.testing.assert_allclose(parray.py(), gdfmap["By"] * momentum / p0)
    np.testing.assert_allclose(parray.p(), (energy - reference_energy) / p0)

    np.testing.assert_allclose(parray.tau(), tau)

    np.testing.assert_allclose(parray.q_array, gdfmap["nmacro"] * abs(gdfmap["q"]))


def test_touts_to_particle_arrays(example_tout):
    parray = gpt_adaptor.tout_to_particle_array(example_tout)
    compare_gdf_dict_with_parray(example_tout, parray)


def test_screens_to_particle_arrays(example_screen):
    parray = gpt_adaptor.screen_to_particle_array(example_screen)
    compare_gdf_dict_with_parray(example_screen, parray)
