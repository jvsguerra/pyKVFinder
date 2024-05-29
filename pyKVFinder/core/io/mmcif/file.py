# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading mmCIF file format.
"""

__name__ = "pyKVFinder.core.io.mmcif"
__all__ = ["read_mmcif"]

import pathlib
import warnings
from typing import Dict, Optional, Union

import numpy
from MDAnalysis import Universe
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from ..vdw import VDW, _lookup_radii, read_vdw


def read_mmcif(
    fn: Union[str, pathlib.Path], vdw: Optional[Dict[str, Dict[str, float]]] = None
) -> numpy.ndarray:
    """Reads mmCIF file into numpy.ndarrays.

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to mmCIF file.
    vdw : Dict[str, Dict[str, float]], optional
        A dictionary containing radii values, by default None. If None, use
        output of ``read_vdw()``.

    Returns
    -------
    atomic : numpy.ndarray
        A numpy array with atomic data (residue number, chain, residue name,
        atom name, xyz coordinates and radius) for each atom.

    Raises
    ------
    TypeError
        `fn` must be a string or a pathlib.Path.

    Note
    ----
    The van der Waals radii file defines the radius values for each atom by
    residue and when not defined, it uses a generic value based on the atom
    type. The function by default loads the built-in van der Waals radii file:
    `vdw.json`.

    See Also
    --------
    read_vdw
    read_xyz
    get_vertices
    get_vertices_from_file
    detect
    constitutional
    hydropathy

    Example
    -------
    With the vdW radii dictionary loaded with ``read_vdw``, we can read a target PDB file into Numpy array (atomic data):

    >>> import os
    >>> import pyKVFinder
    >>> from pyKVFinder import read_pdb
    >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
    >>> atomic = read_pdb(pdb)
    >>> atomic
    array([['13', 'E', 'GLU', ..., '-15.642', '-14.858', '1.824'],
       ['13', 'E', 'GLU', ..., '-14.62', '-15.897', '1.908'],
       ['13', 'E', 'GLU', ..., '-13.357', '-15.508', '1.908'],
       ...,
       ['350', 'E', 'PHE', ..., '18.878', '-9.885', '1.908'],
       ['350', 'E', 'PHE', ..., '17.624', '-9.558', '1.908'],
       ['350', 'E', 'PHE', ..., '19.234', '-13.442', '1.69']],
      dtype='<U32')

    .. warning::
        The function takes the `built-in dictionary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/vdw.dat>`_ when the ``vdw`` argument is not specified. If you wish to use a custom van der Waals radii file, you must read it with ``read_vdw`` as shown earlier and pass it as ``read_pdb(pdb, vdw=vdw)``.
    """
    # Check arguments
    if type(fn) not in [str, pathlib.Path]:
        raise TypeError("`fn` must be a string or a pathlib.Path.")

    # Define default vdw file
    if vdw is None:
        vdw = read_vdw(VDW)

    # Read mmCIF file into a dictionary
    mmcif = MMCIF2Dict(fn)

    # Number of atoms and residues
    n_atoms = len(mmcif["_atom_site.id"])
    n_residues = len({*mmcif["_atom_site.label_seq_id"]})

    # Get residue numbers
    resnums = list(dict.fromkeys(mmcif["_atom_site.label_seq_id"]))

    # Get residue names
    res2index = {
        residue: index for index, residue in enumerate(mmcif["_atom_site.label_seq_id"])
    }
    resnames = [
        mmcif["_atom_site.label_comp_id"][res2index[residue]] for residue in resnums
    ]

    # Get list of resindices for each atom
    res2index = {residue: index for index, residue in enumerate(resnums)}
    resindices = [res2index[residue] for residue in mmcif["_atom_site.label_seq_id"]]

    # Create an empty Universe
    u = Universe.empty(
        n_atoms=n_atoms,
        n_residues=n_residues,
        atom_resindex=resindices,
        trajectory=True,
    )

    # Add topology attributes
    u.add_TopologyAttr(
        "resnums",
        resnums,
    )  # residue number
    u.add_TopologyAttr("chainIDs", mmcif["_atom_site.auth_asym_id"])  # chain identifier
    u.add_TopologyAttr("resnames", resnames)  # residue name
    u.add_TopologyAttr("names", mmcif["_atom_site.label_atom_id"])  # atom name
    u.atoms.positions = numpy.c_[
        numpy.asarray(mmcif["_atom_site.Cartn_x"], dtype=float),
        numpy.asarray(mmcif["_atom_site.Cartn_y"], dtype=float),
        numpy.asarray(mmcif["_atom_site.Cartn_z"], dtype=float),
    ]  # xyz coordinates
    u.add_TopologyAttr("elements", mmcif["_atom_site.type_symbol"])  # atom type

    # Vectorize _lookup_radii function and add radii to topology
    u.add_TopologyAttr(
        "radii",
        numpy.vectorize(_lookup_radii)(
            vdw,
            u.atoms.resnames,
            u.atoms.names,
            u.atoms.elements,
        ),
    )

    # Load atomic data
    atomic = numpy.c_[
        u.atoms.resnums,  # atom number (XYZ file does not group atoms in residues)
        u.atoms.chainIDs,  # chain identifier (XYZ file does not group atoms in chains)
        u.atoms.resnames,  # residue name (XYZ file does not group atoms in residues)
        u.atoms.names,  # atom name
        u.atoms.positions,  # xyz coordinates
        u.atoms.radii,  # atom radius
    ]

    return atomic
