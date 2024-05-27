"""
A subpackage for reading and writing structure related data.

Macromolecular structure files (PDB, PDBx/mmCIF, XYZ, etc.)
to load an :class:`~pyKVFinder.core.structure.Structure` object.
"""

__name__ = "pyKVFinder.core.io"

from .mmcif import *
from .pdb import *
from .pqr import *
from .vdw import *
from .xyz import *
