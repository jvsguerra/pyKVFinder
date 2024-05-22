# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading XYZ chemical file format.
"""

__name__ = "pyKVFinder.core.io.xyz"
__all__ = ["read_xyz"]

import pathlib
from typing import Dict, Optional, Union

import numpy

from ..vdw import VDW, read_vdw


def read_xyz(
    fn: Union[str, pathlib.Path], vdw: Optional[Dict[str, Dict[str, float]]] = None
) -> numpy.ndarray:
    """Reads XYZ file into numpy.ndarrays.

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to XYZ file.
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
    The van der Waals radii file defines the radius values for each atom
    by residue and when not defined, it uses a generic value based on the
    atom type. The function by default loads the built-in van der Waals radii
    file: `vdw.json`.

    See Also
    --------
    read_vdw
    read_pdb
    get_vertices
    get_vertices_from_file
    detect
    constitutional
    hydropathy

    Example
    -------
    With the vdW radii dictionary loaded with ``pyKVFinder.read_vdw``, we can
    read a target XYZ file into Numpy arrays (atomic information and atomic
    coordinates):

    >>> import os
    >>> import pyKVFinder
    >>> from pyKVFinder import read_xyz
    >>> xyz = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.xyz')
    >>> atomic = read_xyz(xyz)
    >>> atominfo
    array([['1', 'A', 'UNK', ..., '-15.642', '-14.858', '1.97'],
       ['2', 'A', 'UNK', ..., '-14.62', '-15.897', '1.66'],
       ['3', 'A', 'UNK', ..., '-13.357', '-15.508', '1.66'],
       ...,
       ['2790', 'A', 'UNK', ..., '18.878', '-9.885', '1.66'],
       ['2791', 'A', 'UNK', ..., '17.624001', '-9.558', '1.66'],
       ['2792', 'A', 'UNK', ..., '19.233999', '-13.442', '1.69']],
      dtype='<U32')

    .. warning::
        The function takes the `built-in dictionary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/core/io/vdw/vdw.json>`_ when the ``vdw`` argument is not specified. If you wish to use a custom van der Waals radii file, you must read it with ``read_vdw`` as shown earlier and pass it as ``read_xyz(xyz, vdw=vdw)``.

    """
    # Check arguments
    if type(fn) not in [str, pathlib.Path]:
        raise TypeError("`fn` must be a string or a pathlib.Path.")

    # Define default vdw file
    if vdw is None:
        vdw = read_vdw(VDW)

    # Create lists
    atomic = []

    # Start resnum
    resnum = 0

    # Read XYZ file
    with open(fn, "r", encoding="utf-8") as f:
        for line in f.readlines():
            line = line.split()
            if len(line) == 4:
                # Get PDB information
                atom_symbol = line[0]
                x = float(line[1])
                y = float(line[2])
                z = float(line[3])

                # Get radius (generic value)
                radius = vdw["GEN"][atom_symbol]

                # Get resnum
                resnum += 1

                # Append data
                atomic.append([resnum, "A", "UNK", atom_symbol, x, y, z, radius])

    return numpy.asarray(atomic)
