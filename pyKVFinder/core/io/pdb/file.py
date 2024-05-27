# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading PDB file format.
"""

__name__ = "pyKVFinder.core.io.pdb"
__all__ = ["read_pdb"]

import pathlib
import warnings
from typing import Dict, Optional, Union

import numpy
from MDAnalysis import Universe

from ..vdw import VDW, _lookup_radii, read_vdw


def read_pdb(
    fn: Union[str, pathlib.Path],
    vdw: Optional[Dict[str, Dict[str, float]]] = None,
    model: Optional[int] = None,
) -> numpy.ndarray:
    """Reads PDB file into numpy.ndarrays.

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to PDB file.
    vdw : Dict[str, Dict[str, float]], optional
        A dictionary containing radii values, by default None. If None, use
        output of ``read_vdw()``.
    model : int, optional
        The model number of a multi-model PDB file, by default None. If None,
        keep atoms from all models.

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
    `vdw.dat`.

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
    if model is not None:
        if type(model) not in [int]:
            raise TypeError("`model` must be an integer.")
        if model < 0:
            raise ValueError("`model` must be a positive integer.")

    # Define default vdw file
    if vdw is None:
        vdw = read_vdw(VDW)

    # Read PDB file in MDAnalysis.Universe
    u = Universe(fn, in_memory=True)
    if model is None:
        if u.trajectory.n_frames > 1:
            warnings.warn(
                f"The {fn} is a multi-model PDB file. Please specify a model \
number (n) in read_pdb({fn}, model=n). Otherwise, the function will read the \
first model."
            )
    else:
        # Select frame number to get positions
        u.trajectory[model]

    # Check if residues names are available
    if "resnames" not in u.atoms._SETATTR_WHITELIST:
        u.add_TopologyAttr("resnames", ["UNK"] * len(u.residues))

    # Check if chain identifiers are available
    if "chainIDs" not in u.atoms._SETATTR_WHITELIST:
        u.add_TopologyAttr("chainIDs", [""])

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
        u.atoms.ids,  # atom number (XYZ file does not group atoms in residues)
        u.atoms.chainIDs,  # chain identifier (XYZ file does not group atoms in chains)
        u.atoms.resnames,  # residue name (XYZ file does not group atoms in residues)
        u.atoms.names,  # atom name
        u.atoms.positions,  # xyz coordinates
        u.atoms.radii,  # atom radius
    ]

    return atomic
