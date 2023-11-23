"""
This module contains methods for reading/writing in openpmd format.
Adapted for Ocelot from APTools:
https://github.com/AngelFP/APtools/blob/master/aptools/data_handling/saving.py
@author: A. Ferran Pousa, A. Martinez de la Ossa.
"""

import os
import numpy as np
import scipy.constants as ct
import h5py

try:
    from openpmd_viewer import OpenPMDTimeSeries
    openpmd_viewer_installed = True
except ImportError:
    openpmd_viewer_installed = False

try:
    from openpmd_api import (Series, Access, Dataset, Mesh_Record_Component,
                             Unit_Dimension)
    openpmd_api_installed = True
except ImportError:
    openpmd_api_installed = False


from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.physics_proc import PhysProc
from ocelot.common.globals import m_e_GeV
from ocelot import __version__

from ocelot.common.ocelog import logging
_logger = logging.getLogger(__name__)


def read_openpmd_to_parray(openpmd_data_path, iteration=None, species=None,
                           select=None, gamma_ref=None, z_ref=None, **kw):
    """
    Reads openpmd particle data.

    Parameters:
    -----------
    openpmd_data_path : Path to openpmd data (dir or file)

    iteration : int
        Iteration number from where the beam data will be read.
        If not specified, it takes the last iteration present in the data.

    gamma_ref : float (Optional)
        Reference energy of the particle beam used for tracking in Ocelot. If
        not specified, the reference energy will be taken as the average
        energy of the input distribution.

    z_ref : float (Optional)
        Reference longitudinal position of the particle beam used for tracking
        in Ocelot. If not specified, the reference value will be taken as the
        average longitudinal position of the input distribution.

    kw: optional arguments to pass to `OpenPMDTimeSeries.get_particle`
        e.g. species, select, ...

    Returns:
    --------
    An Ocelot ParticleArray.

    """

    x, y, z, px, py, pz, q = read_openpmd_to_beam_data(openpmd_data_path,
                                                       iteration=iteration,
                                                       species=species,
                                                       select=select, **kw)
    p_array = beam_data_to_parray([x, y, z, px, py, pz, q],
                                  gamma_ref=gamma_ref, z_ref=z_ref)

    return p_array


def write_parray_to_openpmd(p_array, folder_path, file_name=None,
                            species='beam', iteration=0):
    """
    Writes p_array to an HDF5 file following the openPMD standard.

    Parameters
    ----------
    p_array : ParticleArray
        The Ocelot distribution

    folder_path : str
        Path to the folder in which to save the data

    file_name : str
        Name of the file to save without extension

    species : str
        Optional. Name under which the particle species should be stored.

    iteration : int
        Optional. Iteration number.

    """

    # Get beam data
    beam_data = parray_to_beam_data(p_array)
    write_beam_data_to_openpmd(beam_data, folder_path, file_name=file_name,
                               species=species, iteration=iteration)


def read_openpmd_to_beam_data(openpmd_data_path, iteration=None, species=None,
                              select=None, **kw):
    """
    Reads openpmd particle data.

    Parameters:
    -----------
    openpmd_data_path : Path to openpmd data (dir or file)

    iteration : int
        Iteration number from where the beam data will be read.
        If not specified, it takes the last iteration present in the data.

    species : str
        Optional. Name under which the particle species is stored.

    select : dict
        A dictionary of rules to filter the dataframe, e.g.
        'pz' : [1000, None] (get particles with pz higher than 1000)

    kw: optional arguments to pass to `OpenPMDTimeSeries.get_particle`
        e.g. species, select, ...

    Returns:
    --------
    A list with the beam data as [x, y, z, px, py, pz, q], where the positions
    have units of meters, momentun is in non-dimensional units (beta*gamma)
    and q is in Coulomb.

    """

    if os.path.isdir(openpmd_data_path):
        if not openpmd_viewer_installed:
            raise ImportError('Openpmd_viewer is not installed. '
                              'Cannot read openpmd data folder.')

        diag = OpenPMDTimeSeries(openpmd_data_path,
                                 check_all_files=False, backend='h5py')

        varlist = ['x', 'y', 'z', 'ux', 'uy', 'uz', 'w', 'charge']
        if iteration is None:
            iteration = diag.iterations[-1]
        x, y, z, px, py, pz, w, charge = diag.get_particle(varlist,
                                                           iteration=iteration,
                                                           species=species,
                                                           select=select, **kw)
        q = w * np.fabs(charge)

        return [x, y, z, px, py, pz, q]

    elif os.path.isfile(openpmd_data_path):
        file_content = h5py.File(openpmd_data_path, mode='r')
        # get base path in file
        iteration = list(file_content['/data'].keys())[0]
        base_path = '/data/{}'.format(iteration)
        # get path under which particle data is stored
        particles_path = file_content.attrs['particlesPath'].decode()
        # list available species
        available_species = list(
            file_content[os.path.join(base_path, particles_path)])
        if species is None:
            if len(available_species) == 1:
                species = available_species[0]
            else:
                raise ValueError(
                    'More than one particle species is available. '
                    'Please specify a `species`. '
                    'Available species are: ' + str(available_species))
        # get species
        beam_species = file_content[
            os.path.join(base_path, particles_path, species)]
        # get data
        mass = beam_species['mass']
        charge = beam_species['charge']
        position = beam_species['position']
        position_off = beam_species['positionOffset']
        momentum = beam_species['momentum']
        m = mass.attrs['value'] * mass.attrs['unitSI']
        q = charge.attrs['value'] * charge.attrs['unitSI']
        x = position['x'][:] * position['x'].attrs['unitSI'] \
            + position_off['x'].attrs['value'] \
            * position_off['x'].attrs['unitSI']
        y = position['y'][:] * position['y'].attrs['unitSI'] \
            + position_off['y'].attrs['value'] \
            * position_off['y'].attrs['unitSI']
        z = position['z'][:] * position['z'].attrs['unitSI'] \
            + position_off['z'].attrs['value'] \
            * position_off['z'].attrs['unitSI']
        px = momentum['x'][:] * momentum['x'].attrs['unitSI'] / (m * ct.c)
        py = momentum['y'][:] * momentum['y'].attrs['unitSI'] / (m * ct.c)
        pz = momentum['z'][:] * momentum['z'].attrs['unitSI'] / (m * ct.c)
        w = beam_species['weighting'][:]
        q = w * np.fabs(q)

        if select is not None:
            x, y, z, px, py, pz, q = select_particles([x, y, z, px, py, pz, q],
                                                      select=select)
            print(q)
        return [x, y, z, px, py, pz, q]
    else:
        raise OSError('OpenPMD data not found in ' + str(openpmd_data_path))


def beam_data_to_parray(beam_data, gamma_ref=None, z_ref=None):
    """
    Converts standard beam data to Ocelot ParticleArray

    Parameters:
    -----------
    beam_data : list
        Contains the beam data as [x, y, z, px, py, pz, q], where the positions
        have units of meters, momentun is in non-dimensional units (beta*gamma)
        and q is in Coulomb.

    gamma_ref : float (Optional)
        Reference energy of the particle beam used for tracking in Ocelot. If
        not specified, the reference energy will be taken as the average
        energy of the input distribution.

    z_ref : float (Optional)
        Reference longitudinal position of the particle beam used for tracking
        in Ocelot. If not specified, the reference value will be taken as the
        average longitudinal position of the input distribution.

    Returns:
    --------
    An Ocelot ParticleArray.

    """
    x, y, z, px, py, pz, q = beam_data

    # Calculate gamma.
    gamma = np.sqrt(1 + px**2 + py**2 + pz**2)

    # Determine reference energy and longitudinal position.
    if gamma_ref is None:
        gamma_ref = np.average(gamma, weights=q)
    if z_ref is None:
        z_ref = np.average(z, weights=q)

    # Calculate momentum deviation (dp) and kinetic momentum (p_kin).
    b_ref = np.sqrt(1 - gamma_ref**(-2))
    dp = (gamma - gamma_ref) / (gamma_ref * b_ref)
    p_kin = np.sqrt(gamma_ref**2 - 1)

    # Create particle array
    p_array = ParticleArray(len(q))
    p_array.rparticles[0] = x
    p_array.rparticles[2] = y
    p_array.rparticles[4] = -(z - z_ref)
    p_array.rparticles[1] = px / p_kin
    p_array.rparticles[3] = py / p_kin
    p_array.rparticles[5] = dp
    p_array.q_array[:] = q
    p_array.s = z_ref
    p_array.E = gamma_ref * m_e_GeV

    return p_array


def parray_to_beam_data(p_array):
    """
    Converts an Ocelot ParticleArray to standard beam data

    Parameters:
    -----------
    p_array : ParticleArray
        The Ocelot distribution to be converted.

    Returns:
    --------
    A list with the beam data as [x, y, z, px, py, pz, q], where the positions
    have units of meters, momentun is in non-dimensional units (beta*gamma)
    and q is in Coulomb.
    """

    # Get beam matrix and reference values.
    beam_matrix = p_array.rparticles
    gamma_ref = p_array.E / m_e_GeV
    z_ref = p_array.s

    # Calculate gamma and kinetic momentum (p_kin).
    dp = beam_matrix[5]
    b_ref = np.sqrt(1 - gamma_ref**(-2))
    gamma = dp * gamma_ref * b_ref + gamma_ref
    p_kin = np.sqrt(gamma_ref**2 - 1)

    # Create coordinate arrays
    x = beam_matrix[0]
    px = beam_matrix[1] * p_kin
    y = beam_matrix[2]
    py = beam_matrix[3] * p_kin
    z = z_ref - beam_matrix[4]
    pz = np.sqrt(gamma**2 - px**2 - py**2 - 1)
    q = p_array.q_array

    return [x, y, z, px, py, pz, q]


def write_beam_data_to_openpmd(beam_data, folder_path, file_name=None,
                               species='beam', iteration=0):
    """
    Writes p_array to an HDF5 file following the openPMD standard.

    Parameters
    ----------
        beam_data : list
        Contains the beam data as [x, y, z, px, py, pz, q], where the positions
        have units of meters, momentun is in non-dimensional units (beta*gamma)
        and q is in Coulomb.

    folder_path : str
        Path to the folder in which to save the data

    file_name : str
        Name of the file to save without extension

    species : str
        Optional. Name under which the particle species should be stored.

    iteration : int
        Optional. Iteration number.

    """

    if not openpmd_api_installed:
        raise ImportError('Openpmd_api is not installed. '
                          'Cannot write openpmd data.')

    # Get beam data
    x = np.ascontiguousarray(beam_data[0])
    y = np.ascontiguousarray(beam_data[1])
    z = np.ascontiguousarray(beam_data[2])
    px = np.ascontiguousarray(beam_data[3])
    py = np.ascontiguousarray(beam_data[4])
    pz = np.ascontiguousarray(beam_data[5])
    q = np.ascontiguousarray(beam_data[6])

    # File name
    if file_name is None:
        file_name = 'beam_data%08T.h5'

    # Save to file
    file_path = os.path.join(folder_path, file_name)
    if not file_path.endswith('.h5'):
        file_path += '.h5'

    # OpenPMD series object
    opmd_series = Series(file_path, Access.create)

    # Set basic attributes.
    opmd_series.set_software('ocelot', __version__)
    opmd_series.set_particles_path('particles')

    # Create iteration
    it = opmd_series.iterations[iteration]
    it.time = np.average(z, weights=q) / ct.c

    # Create particles species.
    particles = it.particles[species]

    # Create additional necessary arrays and constants.
    w = np.abs(q) / ct.e
    m = ct.m_e
    q = -ct.e
    px = px * m * ct.c
    py = py * m * ct.c
    pz = pz * m * ct.c

    SCALAR = Mesh_Record_Component.SCALAR

    # Generate datasets.
    d_x = Dataset(x.dtype, extent=x.shape)
    d_y = Dataset(y.dtype, extent=y.shape)
    d_z = Dataset(z.dtype, extent=z.shape)
    d_px = Dataset(px.dtype, extent=px.shape)
    d_py = Dataset(py.dtype, extent=py.shape)
    d_pz = Dataset(pz.dtype, extent=pz.shape)
    d_w = Dataset(w.dtype, extent=w.shape)
    d_q = Dataset(np.dtype('float64'), extent=[1])
    d_m = Dataset(np.dtype('float64'), extent=[1])
    d_xoff = Dataset(np.dtype('float64'), extent=[1])
    d_yoff = Dataset(np.dtype('float64'), extent=[1])
    d_zoff = Dataset(np.dtype('float64'), extent=[1])

    # Record data.
    particles['position']['x'].reset_dataset(d_x)
    particles['position']['y'].reset_dataset(d_y)
    particles['position']['z'].reset_dataset(d_z)
    particles['positionOffset']['x'].reset_dataset(d_xoff)
    particles['positionOffset']['y'].reset_dataset(d_yoff)
    particles['positionOffset']['z'].reset_dataset(d_zoff)
    particles['momentum']['x'].reset_dataset(d_px)
    particles['momentum']['y'].reset_dataset(d_py)
    particles['momentum']['z'].reset_dataset(d_pz)
    particles['weighting'][SCALAR].reset_dataset(d_w)
    particles['charge'][SCALAR].reset_dataset(d_q)
    particles['mass'][SCALAR].reset_dataset(d_m)

    # Prepare for writting.
    particles['position']['x'].store_chunk(x)
    particles['position']['y'].store_chunk(y)
    particles['position']['z'].store_chunk(z)
    particles['positionOffset']['x'].make_constant(0.)
    particles['positionOffset']['y'].make_constant(0.)
    particles['positionOffset']['z'].make_constant(0.)
    particles['momentum']['x'].store_chunk(px)
    particles['momentum']['y'].store_chunk(py)
    particles['momentum']['z'].store_chunk(pz)
    particles['weighting'][SCALAR].store_chunk(w)
    particles['charge'][SCALAR].make_constant(q)
    particles['mass'][SCALAR].make_constant(m)

    # Set units.
    particles['position'].unit_dimension = {Unit_Dimension.L: 1}
    particles['positionOffset'].unit_dimension = {Unit_Dimension.L: 1}
    particles['momentum'].unit_dimension = {
        Unit_Dimension.L: 1,
        Unit_Dimension.M: 1,
        Unit_Dimension.T: -1,
    }
    particles['charge'].unit_dimension = {
        Unit_Dimension.T: 1,
        Unit_Dimension.I: 1,
    }
    particles['mass'].unit_dimension = {Unit_Dimension.M: 1}

    # Set weighting attributes.
    particles['position'].set_attribute('macroWeighted', np.uint32(0))
    particles['positionOffset'].set_attribute(
        'macroWeighted', np.uint32(0))
    particles['momentum'].set_attribute('macroWeighted', np.uint32(0))
    particles['weighting'][SCALAR].set_attribute(
        'macroWeighted', np.uint32(1))
    particles['charge'][SCALAR].set_attribute(
        'macroWeighted', np.uint32(0))
    particles['mass'][SCALAR].set_attribute('macroWeighted', np.uint32(0))
    particles['position'].set_attribute('weightingPower', 0.)
    particles['positionOffset'].set_attribute('weightingPower', 0.)
    particles['momentum'].set_attribute('weightingPower', 1.)
    particles['weighting'][SCALAR].set_attribute('weightingPower', 1.)
    particles['charge'][SCALAR].set_attribute('weightingPower', 1.)
    particles['mass'][SCALAR].set_attribute('weightingPower', 1.)

    # Flush data.
    opmd_series.flush()


def select_particles(beam_data, select=None):
    """
    Filter beam data arrays according to a selection criterium.

    Parameters
    ----------
    beam_data : list
        Contains the beam data as [x, y, z, px, py, pz, w].

    select : dict
        A dictionary of rules to filter the dataframe, e.g.
        'pz' : [1000, None] (get particles with pz higher than 1000)

    Returns
    -------
    A list with the new beam data
    """

    index = {'x': 0, 'y': 1, 'z': 2, 'px': 3, 'py': 4, 'pz': 5, 'w': 6}

    # Create the array that determines whether the particle
    # should be selected or not.
    Ntot = len(beam_data[0])
    select_array = np.ones(Ntot, dtype='bool')

    # Loop over the selection rules,
    # and aggregate results in select_array
    for quantity in select.keys():
        if select[quantity][0] is not None:
            select_array = np.logical_and(
                select_array,
                beam_data[index[quantity]] > select[quantity][0])
        # Check upper bound
        if select[quantity][1] is not None:
            select_array = np.logical_and(
                select_array,
                beam_data[index[quantity]] < select[quantity][1])

    # beam_data = beam_data[:, select_array]
    beam_data_select = []
    for x in beam_data:
        beam_data_select.append(x[select_array])

    return beam_data_select


class SaveBeamOPMD(PhysProc):
    def __init__(self, folder_path, file_name=None,
                 species='beam',
                 iteration=0, **kw):
        # PhysProc.__init__(self)
        super(SaveBeamOPMD, self).__init__(**kw)
        self.folder_path = folder_path
        self.file_name = file_name
        self.species = species
        self.iteration = iteration

    def apply(self, p_array, dz):
        _logger.debug(" SaveBeamOPMD applied")

        write_parray_to_openpmd(p_array=p_array,
                                folder_path=self.folder_path,
                                file_name=self.file_name,
                                species=self.species,
                                iteration=self.iteration)

        self.iteration = self.iteration + 1
