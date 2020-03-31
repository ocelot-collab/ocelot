__author__ = 'Sergey Tomin'

import os
import multiprocessing

num_thread_default = str(int(multiprocessing.cpu_count()/2))


OCELOT_NUM_THREADS = os.environ.get("OCELOT_NUM_THREADS", num_thread_default)
os.environ["NUMEXPR_NUM_THREADS"] = OCELOT_NUM_THREADS
os.environ["NUMBA_NUM_THREADS"] = OCELOT_NUM_THREADS
