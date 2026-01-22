"""Python-C parallel KVFinder.

pyKVFinder detects and characterizes cavities in biomolecular structures.

The characterization includes shape, volume, area, depth, hydropathy and
interface residues and their frequencies.

In addition to the set of functions that can be imported into Python scripts,
it contains a command line interface (CLI).

Python package
--------------
>>> import pyKVFinder

Command Line Interface
----------------------
  Usage: pyKVFinder [-h] [-v] [--version] [-b <str>] [-O <str>]\
                    [-m <int>] [--nthreads <int>] [-d <str>] [-s <float>]\
                    [-i <float>] [-o <float>] [-V <float>] [-R <float>]\
                    [-S <str>] [--ignore_backbone] [-D] [--plot_frequencies]\
                    [--hydropathy [{EisenbergWeiss, HessaHeijne, KyteDoolittle,\
                    MoonFleming, RadzickaWolfenden, WimleyWhite, ZhaoLondon, <.toml>}]]\
                    [-B <.toml>] [-L (<.pdb> | <.xyz>)] [--ligand_cutoff <float>]\
                    (<.pdb> | <.xyz>)

See also
--------
* GitHub repository: https://github.com/LBC-LNBio/pyKVFinder

* Documentation: https://lbc-lnbio.github.io/pyKVFinder
"""

__name__ = "pyKVFinder"
__version__ = "0.9.0"
__license__ = "GNU GPL-3.0 License"

from .grid import (
    constitutional,
    depth,
    detect,
    export,
    export_openings,
    get_vertices,
    get_vertices_from_file,
    hydropathy,
    openings,
    spatial,
)
from .main import Molecule, pyKVFinderResults, run_workflow
from .utils import (
    calculate_frequencies,
    plot_frequencies,
    read_cavity,
    read_pdb,
    read_vdw,
    read_xyz,
    write_results,
)
