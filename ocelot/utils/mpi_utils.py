__author__ = 'Sergey Tomin'

from ocelot.cpbd.beam import ParticleArray
import numpy as np


def slicing_p_array(p_array, SIZE):
    """
    Divide ParticleArray into chuncks

    :param p_array: ParticleArray
    :param SIZE: number of threads
    :return: list of ParticleArray
    """

    if SIZE == 0:
        raise ValueError("SIZE must be not zero")
    if SIZE > p_array.n:
        raise ValueError("SIZE must be less or equal number of paricles")
    p_array_list = []
    nchunk = int(p_array.n / SIZE)
    nparticles_list = np.ones(SIZE, dtype=int) * nchunk
    nparticles_last = p_array.n - nchunk * (SIZE - 1)
    if nparticles_last > nchunk:
        nchunk += 1
        nextra = nparticles_last - nchunk
        nparticles_list[:nextra] += 1
        nparticles_list[-1] = nparticles_last - nextra
    ncontrol = 0

    for i, nparticles in enumerate(nparticles_list):
        p_array_part = ParticleArray(n=nparticles)
        p_array_part.rparticles[:, :] = p_array.rparticles[:, ncontrol:ncontrol + nparticles]
        p_array_part.q_array[:] = p_array.q_array[ncontrol:ncontrol + nparticles]
        ncontrol += nparticles
        p_array_part.s = p_array.s
        p_array_part.E = p_array.E
        p_array_list.append(p_array_part)
    return p_array_list


def chunks(lst, nthreads):
    """
    Divide list into chuncks

    :param lst: list of something
    :param nthreads: number of threads
    :return: list of list
    """
    if len(lst) == 0:
        raise ValueError("Length of list must be not zero")
    if nthreads == 0:
        raise ValueError("Number of threads must be not zero")
    if nthreads > len(lst):
        raise ValueError("Number of threads ({}) must be less or equal list length ({})".format(nthreads, len(lst)))

    n = int(len(lst) / nthreads)
    return [lst[i:i + n] for i in range(0, len(lst), n)]


if __name__ == "__main__":
    p_array = ParticleArray(n=10)
    p_array_list = slicing_p_array(p_array, SIZE=1)
    print("test1: len(list) = ", len(p_array_list), "n particles in each: ", [pa.n for pa in p_array_list])

    p_array_list = slicing_p_array(p_array, SIZE=2)
    print("test2: len(list) = ", len(p_array_list), "n particles in each: ", [pa.n for pa in p_array_list])

    p_array_list = slicing_p_array(p_array, SIZE=3)
    print("test3: len(list) = ", len(p_array_list), "n particles in each: ", [pa.n for pa in p_array_list])

    p_array_list = slicing_p_array(p_array, SIZE=10)
    print("test4: len(list) = ", len(p_array_list), "n particles in each: ", [pa.n for pa in p_array_list])

    p_array_list = slicing_p_array(p_array, SIZE=0)
    print("test6: len(list) = ", len(p_array_list), "n particles in each: ", [pa.n for pa in p_array_list])

    p_array_list = slicing_p_array(p_array, SIZE=11)
    print("test5: len(list) = ", len(p_array_list), "n particles in each: ", [pa.n for pa in p_array_list])
