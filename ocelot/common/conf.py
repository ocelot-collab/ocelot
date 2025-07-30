__author__ = 'Sergey Tomin'

import os
import multiprocessing

num_thread_default = str(int(multiprocessing.cpu_count()/2))

#num_thread_default = "1"
OCELOT_NUM_THREADS = os.environ.get("OCELOT_NUM_THREADS", num_thread_default)
os.environ["NUMEXPR_NUM_THREADS"] = OCELOT_NUM_THREADS
os.environ["NUMBA_NUM_THREADS"] = OCELOT_NUM_THREADS
#os.environ["OMP_NUM_THREADS"] = OCELOT_NUM_THREADS # export OMP_NUM_THREADS=1
#os.environ["MKL_NUM_THREADS"] = OCELOT_NUM_THREADS # export MKL_NUM_THREADS=1

"""
export OCELOT_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
"""