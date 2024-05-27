# This source code is part of the pyKVFinder package and is distributed
# under the GNU GPL-3.0 license. Please see 'LICENSE' for further
# information.

"""
This subpackage is used for reading van der Waals radii file.
"""

__name__ = "pyKVFinder.core.io.vdw"
__all__ = ["VDW", "read_vdw"]

import json
import os
import pathlib
from typing import Dict, Optional, Union

VDW = os.path.join(os.path.abspath(os.path.dirname(__file__)), "vdw.json")


def _lookup_radii(vdw: Dict[str, Dict[str, float]], resname: str, atomname: str, atomsymbol: str) -> float:
    """
    Look up the van der Waals radius for a given atom in a specific residue.

    Parameters
    ----------
    vdw : Dict[str, Dict[str, float]]
        A dictionary containing radii values.
    resname : str
        The name of the residue.
    atomname : str
        The name of the atom.
    atomsymbol : str
        The element symbol of the atom.

    Returns
    -------
    float
        The van der Waals radius of the atom. If the specific radius is not
        found, the function falls back to generic values.

    Raises
    ------
    KeyError
        If the radius for the atom or symbol is not found in the provided dictionary.
        - If `resname` exists in `vdw` but `atomname` is not found in the corresponding dictionary.
        - If `resname` does not exist in `vdw` and `atomsymbol` is not found in the 'GEN' category.
    """
    try:
        # Look up the radius for the specific residue and atom
        return vdw.get(resname, vdw["GEN"]).get(
            atomname.upper(), vdw["GEN"][atomsymbol.upper()]
        )
    except KeyError as e:
        if resname in vdw and atomname.upper() not in vdw[resname]:
            raise KeyError(
                f"Van der Waals radius for atom '{atomname}' not found in residue '{resname}'"
            ) from e
        if resname not in vdw and atomsymbol.upper() not in vdw["GEN"]:
            raise KeyError(
                f"Van der Waals radius for atom symbol '{atomsymbol}' not found in 'GEN' category"
            ) from e
        raise e


def read_vdw(
    fn: Optional[Union[str, pathlib.Path]] = None
) -> Dict[str, Dict[str, float]]:
    """Reads van der Waals radii from .json file.

    Parameters
    ---------->
    fn : Optional[Union[str, pathlib.Path]], optional
        A path to a van der Waals radii file, by default None. If None, apply the built-in van der
        Waals radii file `vdw.json`.

    Returns
    -------
    vdw : Dict[str, Dict[str, float]]
        A dictionary containing radii values.

    Raises
    ------
    TypeError
        `fn` must be a string or a pathlib.Path.

    Note
    ----
    The van der Waals radii file defines the radius values for each
    atom by residue and when not defined, it uses a generic value
    based on the atom type (see ``van der Waals file template``).
    The package contains a built-in van der Waals radii file: ``vdw.json``.

    See Also
    --------
    read_pdb
    read_xyz
    Molecule

    Example
    -------
    The ``read_vdw`` function takes the `built-in dictionary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/core/io/vdw/vdw.json>`_ when a *.json* file is not specified.

    >>> from pyKVFinder import read_vdw
    >>> vdw = read_vdw()
    >>> vdw
    {'ALA': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB1': 1.487, '1HB': 1.487, 'HB2': 1.487, '2HB': 1.487, 'HB3': 1.487, '3HB': 1.487, 'C': 1.908, 'O': 1.6612}, 'ARG': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'HD2': 1.387, '1HD': 1.387, '2HD': 1.387, 'HD3': 1.387, 'HD1': 1.387, 'NE': 1.75, 'HE': 0.6, 'CZ': 1.908, 'NH1': 1.75, 'HH11': 0.6, '1HH1': 0.6, 'HH12': 0.6, '2HH1': 0.6, 'NH2': 1.75, 'HH21': 0.6, '2HH2': 0.6, 'HH22': 0.6, '1HH2': 0.6, 'C': 1.908, 'O': 1.6612}, 'ASH': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'OD2': 1.721, 'HD2': 0.0001, 'C': 1.908, 'O': 1.6612}, 'ASN': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'ND2': 1.824, 'HD21': 0.6, '1HD2': 0.6, 'HD22': 0.6, '2HD2': 0.6, 'C': 1.908, 'O': 1.6612}, 'ASP': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'OD2': 1.6612, 'C': 1.908, 'O': 1.6612}, 'CYM': {'N': 1.824, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB3': 1.387, 'HB2': 1.387, 'SG': 2.0, 'C': 1.908, 'O': 1.6612}, 'CYS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, '2HB': 1.387, '1HB': 1.387, 'HB3': 1.387, 'HB1': 1.387, 'SG': 2.0, 'HG': 0.6, 'C': 1.908, 'O': 1.6612}, 'CYX': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, 'SG': 2.0, 'C': 1.908, 'O': 1.6612}, 'GLH': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'HG2': 1.487, 'HG3': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'OE2': 1.721, 'HE2': 0.0001, 'C': 1.908, 'O': 1.6612}, 'GLN': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'NE2': 1.824, 'HE21': 0.6, '1HE2': 0.6, 'HE22': 0.6, '2HE2': 0.6, 'C': 1.908, 'O': 1.6612}, 'GLU': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'OE2': 1.6612, 'C': 1.908, 'O': 1.6612}, 'GLY': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA2': 1.387, 'HA1': 1.387, '1HA': 1.387, '2HA': 1.387, 'HA3': 1.387, 'C': 1.908, 'O': 1.6612}, 'HID': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'HIE': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'HIP': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'ILE': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.487, 'CG2': 1.908, 'HG21': 1.487, '1HG2': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG23': 1.487, '3HG2': 1.487, 'CG1': 1.908, 'HG12': 1.487, '2HG1': 1.487, 'HG13': 1.487, 'HG11': 1.487, '1HG1': 1.487, 'CD1': 1.908, 'HD11': 1.487, '1HD1': 1.487, 'HD12': 1.487, '2HD1': 1.487, 'HD13': 1.487, '3HD1': 1.487, 'C': 1.908, 'O': 1.6612}, 'LEU': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG': 1.487, 'CD1': 1.908, 'HD11': 1.487, '1HD1': 1.487, 'HD12': 1.487, '2HD1': 1.487, 'HD13': 1.487, '3HD1': 1.487, 'CD2': 1.908, 'HD21': 1.487, '1HD2': 1.487, 'HD22': 1.487, '2HD2': 1.487, 'HD23': 1.487, '3HD2': 1.487, 'C': 1.908, 'O': 1.6612}, 'LYN': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'HG2': 1.487, 'HG3': 1.487, 'CD': 1.908, 'HD2': 1.487, 'HD3': 1.487, 'CE': 1.908, 'HE2': 1.1, 'HE3': 1.1, 'NZ': 1.824, 'HZ2': 0.6, 'HZ3': 0.6, 'C': 1.908, 'O': 1.6612}, 'LYS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'HD2': 1.487, '1HD': 1.487, '2HD': 1.487, 'HD3': 1.487, 'HD1': 1.487, 'CE': 1.908, 'HE2': 1.1, '2HE': 1.1, 'HE3': 1.1, '1HE': 1.1, 'HE1': 1.1, 'NZ': 1.824, 'HZ1': 0.6, '1HZ': 0.6, 'HZ2': 0.6, '2HZ': 0.6, 'HZ3': 0.6, '3HZ': 0.6, 'C': 1.908, 'O': 1.6612}, 'MET': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.387, '2HG': 1.387, 'HG3': 1.387, 'HG1': 1.387, '1HG': 1.387, 'SD': 2.0, 'CE': 1.908, 'HE1': 1.387, '1HE': 1.387, 'HE2': 1.387, '2HE': 1.387, 'HE3': 1.387, '3HE': 1.387, 'C': 1.908, 'O': 1.6612}, 'PHE': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'HZ': 1.459, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'C': 1.908, 'O': 1.6612}, 'PRO': {'N': 1.824, 'CD': 1.908, 'HD2': 1.387, '1HD': 1.387, '2HD': 1.387, 'HD3': 1.387, 'HD1': 1.387, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CA': 1.908, 'HA': 1.387, 'C': 1.908, 'O': 1.6612}, 'SER': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, '2HB': 1.387, '1HB': 1.387, 'HB3': 1.387, 'HB1': 1.387, 'OG': 1.721, 'HG': 0.0001, 'C': 1.908, 'O': 1.6612}, 'THR': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, '1HG2': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG23': 1.487, '3HG2': 1.487, 'OG1': 1.721, 'HG1': 0.0001, 'C': 1.908, 'O': 1.6612}, 'TRP': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.85, 'CD1': 2.0, 'HD1': 1.409, 'NE1': 1.75, 'HE1': 0.6, 'CE2': 1.85, 'CZ2': 1.908, 'HZ2': 1.459, 'CH2': 1.908, 'HH2': 1.459, 'CZ3': 1.908, 'HZ3': 1.459, 'CE3': 1.908, 'HE3': 1.459, 'CD2': 1.85, 'C': 1.908, 'O': 1.6612}, 'TYR': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'OH': 1.721, 'HH': 0.0001, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'C': 1.908, 'O': 1.6612}, 'VAL': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.487, 'CG1': 1.908, 'CG2': 1.908, 'HG11': 1.487, '1HG2': 1.487, '1HG1': 1.487, 'HG21': 1.487, 'HG12': 1.487, '2HG1': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG13': 1.487, '3HG2': 1.487, '3HG1': 1.487, 'HG23': 1.487, 'C': 1.908, 'O': 1.6612}, 'HIS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'PTR': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'OH': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'SEP': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, '1HB': 1.387, '2HB': 1.387, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'TPO': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, 'HG22': 1.487, 'HG23': 1.487, '1HG2': 1.487, '2HG2': 1.487, '3HG2': 1.487, 'OG1': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'H2D': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'Y1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'T1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, 'HG22': 1.487, 'HG23': 1.487, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'S1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'GEN': {'AC': 2.0, 'AG': 1.72, 'AL': 2.0, 'AM': 2.0, 'AR': 1.88, 'AS': 1.85, 'AT': 2.0, 'AU': 1.66, 'B': 2.0, 'BA': 2.0, 'BE': 2.0, 'BH': 2.0, 'BI': 2.0, 'BK': 2.0, 'BR': 1.85, 'C': 1.66, 'CA': 2.0, 'CD': 1.58, 'CE': 2.0, 'CF': 2.0, 'CL': 1.75, 'CM': 2.0, 'CO': 2.0, 'CR': 2.0, 'CS': 2.0, 'CU': 1.4, 'DB': 2.0, 'DS': 2.0, 'DY': 2.0, 'ER': 2.0, 'ES': 2.0, 'EU': 2.0, 'F': 1.47, 'FE': 2.0, 'FM': 2.0, 'FR': 2.0, 'GA': 1.87, 'GD': 2.0, 'GE': 2.0, 'H': 0.91, 'HE': 1.4, 'HF': 2.0, 'HG': 1.55, 'HO': 2.0, 'HS': 2.0, 'I': 1.98, 'IN': 1.93, 'IR': 2.0, 'K': 2.75, 'KR': 2.02, 'LA': 2.0, 'LI': 1.82, 'LR': 2.0, 'LU': 2.0, 'MD': 2.0, 'MG': 1.73, 'MN': 2.0, 'MO': 2.0, 'MT': 2.0, 'N': 1.97, 'NA': 2.27, 'NB': 2.0, 'ND': 2.0, 'NE': 1.54, 'NI': 1.63, 'NO': 2.0, 'NP': 2.0, 'O': 1.69, 'OS': 2.0, 'P': 2.1, 'PA': 2.0, 'PB': 2.02, 'PD': 1.63, 'PM': 2.0, 'PO': 2.0, 'PR': 2.0, 'PT': 1.72, 'PU': 2.0, 'RA': 2.0, 'RB': 2.0, 'RE': 2.0, 'RF': 2.0, 'RH': 2.0, 'RN': 2.0, 'RU': 2.0, 'S': 2.09, 'SB': 2.0, 'SC': 2.0, 'SE': 1.9, 'SG': 2.0, 'SI': 2.1, 'SM': 2.0, 'SN': 2.17, 'SR': 2.0, 'TA': 2.0, 'TB': 2.0, 'TC': 2.0, 'TE': 2.06, 'TH': 2.0, 'TI': 2.0, 'TL': 1.96, 'TM': 2.0, 'U': 1.86, 'V': 2.0, 'W': 2.0, 'XE': 2.16, 'Y': 2.0, 'YB': 2.0, 'ZN': 1.39, 'ZR': 2.0}}

    The van der Waals radii can be define by:

        * creating a Python dictionary:

        >>> vdw = {'GEN': {'C': 1.66, 'CA': 2.0, 'N': 1.97, 'O': 1.69, 'H': 0.91}}

        * specifying a *.json* file following template of `van der Waals radii file`.

        >>> import json
        >>> with open('vdw.json', 'w') as f:
        ...     json.dump(vdw, f)
        >>> vdw = read_vdw('vdw.json')
        >>> vdw
        {'GEN': {'C': 1.66, 'CA': 2.0, 'N': 1.97, 'O': 1.69, 'H': 0.91}}
    """
    # Check argument
    if fn is not None:
        if type(fn) not in [str, pathlib.Path]:
            raise TypeError("`fn` must be a string or a pathlib.Path.")
    else:
        # Define default vdw file
        fn = VDW

    # Check if file exists
    with open(fn, "r", encoding="utf-8") as f:
        vdw = json.load(f)

    return vdw
