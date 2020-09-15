import numpy as np

from ocelot.cpbd.beam import ParticleArray


def test_lost_particle_recorder():
    n = 20
    p_array = np.arange(n)
    lpr = ParticleArray.LostParticleRecorder(n)

    for pos, inds in enumerate([[1, 3], [0, 1], [4, 10], [0, 1], [1, 0, 10, 7, 8]]):
        lpr.add(inds, pos)
        p_array = np.delete(p_array, inds)

    lpr.lost_particles.sort()
    assert all([lp_idx == truth_idx for lp_idx, truth_idx in
                zip(lpr.lost_particles, [0, 1, 2, 3, 4, 5, 6, 7, 8, 14, 15, 16, 18])])

    assert all([lp_idx == truth_idx for lp_idx, truth_idx in
                zip(lpr.lp_to_pos_hist, [(0, 2), (1, 2), (2, 2), (3, 2), (4, 5)])])
