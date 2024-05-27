# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading PQR file format.
"""

__name__ = "pyKVFinder.core.io.pqr"
__all__ = ["read_pqr"]

import pathlib
import warnings
from typing import Union

import numpy
from MDAnalysis import Universe

from ..vdw import VDW, _lookup_radii, read_vdw


def read_pqr(fn: Union[str, pathlib.Path]) -> numpy.ndarray:
    """Reads PQR file into numpy.ndarrays.

    Parameters
    ----------
    fn : Union[str, pathlib.Path]
        A path to PQR file.

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
    read_xyz
    read_pdb
    get_vertices
    get_vertices_from_file
    detect
    constitutional
    hydropathy
    """
    # Check arguments
    if type(fn) not in [str, pathlib.Path]:
        raise TypeError("`fn` must be a string or a pathlib.Path.")

    # Read PQR file in MDAnalysis.Universe
    u = Universe(fn, in_memory=True)
    if u.trajectory.n_frames > 1:
        raise ValueError("The PQR file must contain only one frame.")

    # Check if residues names are available
    if "resnames" not in u.atoms._SETATTR_WHITELIST:
        u.add_TopologyAttr("resnames", ["UNK"] * len(u.residues))

    # Check if chain identifiers are available
    if "chainIDs" not in u.atoms._SETATTR_WHITELIST:
        u.add_TopologyAttr("chainIDs", [""] * len(u.atoms))

    # Check if radii are available
    if "radii" not in u.atoms._SETATTR_WHITELIST:
        # Vectorize _lookup_radii function and add radii to topology
        vdw = read_vdw(VDW)
        u.add_TopologyAttr(
            "radii",
            numpy.vectorize(_lookup_radii)(
                vdw,
                u.atoms.resnames,
                u.atoms.names,
                u.atoms.elements,
            ),
        )
    if u.atoms.radii.sum() == 0.0:
        warnings.warn(
            "The PQR file does not contain radii values. Loading default van \
der Waals radii into PQR data."
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
