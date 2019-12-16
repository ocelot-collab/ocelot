"""
Adaptors to translate CSRtrack particle distribution .fmt1 to ParticleArray and back

S.Tomin and I.Zagorodnov
"""
from ocelot.common.globals import *
from ocelot.cpbd.beam import ParticleArray

def csrtrackBeam2particleArray(filename, orient="H"):
    """
    Function to read CSRtrack beam files ".fmt1"
    H: z x y pz px py -> x y z px py pz
    V: z y x pz py px -> x y -z px py -pz

    :param filename: filename
    :param orient: str, "H" or "V" horizontal or vertical orientation
    :return: ParticleArray
    """
    pd = np.loadtxt(filename)
    n = np.shape(pd)[0] - 1
    pd1 = np.zeros((n, 6))

    if orient == 'H':
        pd1[:, 0] = pd[1:, 1]
        pd1[:, 1] = pd[1:, 2]
        pd1[:, 2] = pd[1:, 0]
        pd1[:, 3] = pd[1:, 4]
        pd1[:, 4] = pd[1:, 5]
        pd1[:, 5] = pd[1:, 3]
    else:
        pd1[:, 0] = -pd[1:, 2]
        pd1[:, 1] = pd[1:, 1]
        pd1[:, 2] = pd[1:, 0]
        pd1[:, 3] = -pd[1:, 5]
        pd1[:, 4] = pd[1:, 4]
        pd1[:, 5] = pd[1:, 3]

    for i in range(6):
        pd1[1:n, i] = pd1[1:n, i] + pd1[0, i]

    p_ref = np.sqrt(pd1[0, 3]**2 + pd1[0, 4]**2 + pd1[0, 5]**2)

    px = pd1[:, 3] / p_ref
    py = pd1[:, 4] / p_ref
    Eref = np.sqrt(m_e_eV ** 2 + p_ref ** 2)
    pe = (np.sqrt(m_e_eV**2 + (pd1[:, 3]**2 + pd1[:, 4]**2 + pd1[:, 5]**2)) - Eref) / p_ref

    p_array = ParticleArray(n)
    p_array.rparticles[0] = pd1[:, 0]
    p_array.rparticles[2] = pd1[:, 1]
    p_array.rparticles[4] = -(pd1[:, 2] - pd1[0, 2])
    p_array.rparticles[1] = px[:]
    p_array.rparticles[3] = py[:]
    p_array.rparticles[5] = pe[:]

    p_array.q_array[:] = pd[1:, 6]
    p_array.s = pd1[0, 2]
    p_array.E = Eref*1e-9
    return p_array

def particleArray2csrtrackBeam(p_array, filename="csr_beam.fmt1"):
    """
    Translate ParticleArray to CSRtrack particle distribution .fmt1

    :param p_array: ParticleArray
    :param filename: filename
    :return:
    """
    Eref = p_array.E * 1e9  # GeV -> eV
    p_ref = np.sqrt(Eref ** 2 - m_e_eV ** 2)
    E = p_array.rparticles[5] * p_ref + Eref
    p = np.sqrt(E**2 - m_e_eV**2)

    px = p_array.rparticles[1] * p_ref
    py = p_array.rparticles[3] * p_ref
    pz = np.sqrt(p ** 2 - px ** 2 - py ** 2)
    x = p_array.rparticles[0]
    y = p_array.rparticles[2]
    z = -p_array.rparticles[4] + p_array.s

    pd = np.zeros((p_array.n + 1, 7))
    t = 0

    pd[0] = [t, 0, 0, 0, 0, 0, 0 ]
    pd[1] = [z[0], x[0], y[0], pz[0], px[0], py[0], p_array.q_array[0]]
    pd[2:, 0] = z[1:] - z[0]
    pd[2:, 1] = x[1:] - x[0]
    pd[2:, 2] = y[1:] - y[0]
    pd[2:, 3] = pz[1:] - pz[0]
    pd[2:, 4] = px[1:] - px[0]
    pd[2:, 5] = py[1:] - py[0]
    pd[2:, 6] = p_array.q_array[1:]
    np.savetxt(filename, pd, fmt='%.7e')