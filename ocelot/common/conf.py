"""Thread-count defaults for optional numerical backends used by Ocelot."""

import multiprocessing
import os
import sys


num_thread_default = str(max(1, multiprocessing.cpu_count() // 2))
OCELOT_NUM_THREADS = os.environ.get("OCELOT_NUM_THREADS", num_thread_default)

# Numba reads its pool size at import. Changing the environment later can
# invalidate a pool already running for another library in the same process.
if "numba" not in sys.modules:
    os.environ.setdefault("NUMBA_NUM_THREADS", OCELOT_NUM_THREADS)

# NumExpr also reads its defaults at import. Leave it to apply the user's
# maximum pool size or OpenMP settings when either has been configured.
if ("numexpr" not in sys.modules
        and not os.environ.get("NUMEXPR_MAX_THREADS")
        and not os.environ.get("OMP_NUM_THREADS")):
    os.environ.setdefault("NUMEXPR_NUM_THREADS", OCELOT_NUM_THREADS)
