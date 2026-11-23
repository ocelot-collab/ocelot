"""Check import order using fresh processes because thread pools are process-wide."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import textwrap

import pytest


def run_with_threads(code: str, **settings: str) -> None:
    env = os.environ.copy()
    for name in ("OCELOT_NUM_THREADS", "NUMBA_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                 "NUMEXPR_MAX_THREADS", "OMP_NUM_THREADS"):
        env.pop(name, None)
    env.update(settings)
    subprocess.run([sys.executable, "-c", textwrap.dedent(code)], env=env,
                   cwd=Path(__file__).resolve().parents[1], check=True)


def test_backend_thread_defaults_and_explicit_overrides():
    run_with_threads("""
        import multiprocessing, os, sys
        from unittest.mock import patch
        with patch.object(multiprocessing, "cpu_count", return_value=1):
            from ocelot.common import conf
        assert conf.OCELOT_NUM_THREADS == "1"
        assert os.environ["NUMBA_NUM_THREADS"] == "1"
        assert os.environ["NUMEXPR_NUM_THREADS"] == "1"
        assert "numba" not in sys.modules and "numexpr" not in sys.modules
    """)
    run_with_threads("""
        import os
        from ocelot.common import conf
        assert conf.OCELOT_NUM_THREADS == "3"
        assert os.environ["NUMBA_NUM_THREADS"] == "2"
        assert os.environ["NUMEXPR_NUM_THREADS"] == "1"
    """, OCELOT_NUM_THREADS="3", NUMBA_NUM_THREADS="2", NUMEXPR_NUM_THREADS="1")


@pytest.mark.parametrize("setting", ["OMP_NUM_THREADS", "NUMEXPR_MAX_THREADS"])
def test_respects_numexpr_thread_settings(setting):
    run_with_threads("""
        import os
        from ocelot.common import conf
        assert "NUMEXPR_NUM_THREADS" not in os.environ
        assert os.environ["NUMBA_NUM_THREADS"] == "3"
    """, OCELOT_NUM_THREADS="3", **{setting: "2"})


@pytest.mark.skipif(importlib.util.find_spec("numba") is None, reason="Numba is optional")
def test_numba_pool_initialized_before_ocelot():
    run_with_threads("""
        import os
        import numba
        threads = numba.get_num_threads()  # Launch the pool with its own defaults.
        os.environ["OCELOT_NUM_THREADS"] = str(1 if threads > 1 else 2)
        from ocelot.common import conf
        assert "NUMBA_NUM_THREADS" not in os.environ
        @numba.njit
        def increment(x):
            return x + 1
        assert increment(2) == 3
        assert numba.get_num_threads() == threads
    """)


@pytest.mark.skipif(importlib.util.find_spec("numexpr") is None, reason="NumExpr is optional")
def test_numexpr_pool_initialized_before_ocelot():
    run_with_threads("""
        import os
        import numexpr
        threads = numexpr.get_num_threads()
        from ocelot.common import conf
        assert "NUMEXPR_NUM_THREADS" not in os.environ
        assert numexpr.get_num_threads() == threads
        assert numexpr.evaluate("1+1") == 2
    """, OCELOT_NUM_THREADS="1")


def test_space_charge_honors_ocelot_threads_when_it_loads_the_backends():
    run_with_threads("""
        import os, sys
        from ocelot.cpbd.sc import SpaceCharge
        assert os.environ["NUMBA_NUM_THREADS"] == "1"
        assert os.environ["NUMEXPR_NUM_THREADS"] == "1"
        for name in ("numba", "numexpr"):
            if name in sys.modules:
                assert sys.modules[name].get_num_threads() == 1
    """, OCELOT_NUM_THREADS="1")
