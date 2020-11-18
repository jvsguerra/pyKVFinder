import os
import sys
import time
import argparse

from . import _name, _version

def positive_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0:
        raise argparse.ArgumentTypeError("%r not a positive floating-point" % (x,))
    return x

def restricted_step_size(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if not 0.0 < x <= 2.0:
        raise argparse.ArgumentTypeError("%r not in range (0.0, 2.0]"%(x,))
    return x

def argparser(__title__: str = None, __version__: str = None, __license__: str = None
    ):

    # Overrides method in HelpFormatter
    class CapitalisedHelpFormatter(argparse.HelpFormatter):

        def add_usage(self, usage, actions, groups, prefix=None):
            if prefix is None:
                prefix = 'Usage: '
            return super(CapitalisedHelpFormatter, self).add_usage(
                usage, actions, groups, prefix)

    parser = argparse.ArgumentParser(
        prog="pyKVFinder",
        description="pyKVFinder software detects and characterizes cavities.",
        formatter_class=CapitalisedHelpFormatter,
        add_help=True)

    # Change parser titles
    parser._positionals.title = 'Positional arguments'
    parser._optionals.title = 'Optional arguments'

    # Positional arguments
    parser.add_argument(
        "pdb",
        metavar="<.pdb>",
        help="Path to a target PDB file."
    )

    # Optional arguments
    parser.add_argument("-v", "--verbose",
                        help="Print extra information to standard output.",
                        action="store_true")
    parser.add_argument('--version',
                        action='version',
                        version=f'{_name} v{_version}',
                        help="Show pyKVFinder version number and exit.")
    parser.add_argument("--base_name",
                        metavar="<str>",
                        help="Prefix to output files.")
    parser.add_argument("--output_directory",
                        metavar="<path>",
                        default=os.path.abspath("."),
                        help="Output directory. (default: %(default)s)")

    # Create argument group
    parameters = parser.add_argument_group("Parameters")

    parameters.add_argument("-d", "--dictionary",
                             default=os.path.join(os.path.abspath(os.path.dirname(__file__)),"data/vdw.dat"), # FIXME: prepare it latter
                             metavar="<file>",
                             type=str,
                             help="Path to van der Waals radii file. (default: %(default)s)")
    parameters.add_argument("-s", "--step",
                             metavar="<float>",
                             default=0.6,
                             type=restricted_step_size,
                             help="Step size (grid spacing). (default: %(default).1lf)")
    parameters.add_argument("-i", "--probe_in",
                             metavar="<float>",
                             default=1.4,
                             type=positive_float,
                             help=u"Probe In size (\u212B). (default: %(default).1lf)")
    parameters.add_argument("-o", "--probe_out",
                             metavar="<float>",
                             default=4.0,
                             type=positive_float,
                             help=u"Probe Out size (\u212B). (default: %(default).1lf)")
    parameters.add_argument("-V", "--volume_cutoff",
                             metavar="<float>",
                             default=5.0,
                             type=positive_float,
                             help=u"Cavities volume filter (\u212B\u00b3). (default: %(default).1lf)")
    parameters.add_argument("-R", "--removal_distance",
                             metavar="<float>",
                             default=2.4,
                             type=positive_float,
                             help=u"Length to be removed from the cavity-bulk frontier (\u212B). \
                             (default: %(default).1lf)")
    parameters.add_argument("-S", "--surface",
                             metavar="<enum>",
                             default="SES",
                             choices=["SES", "SAS"],
                             help="Surface representation. Options: %(choices)s. SAS specifies solvent accessible surface. SES specifies solvent excluded surface. (default: %(default)s)")

    # Create argument group
    box_adjustment = parser.add_argument_group("Box adjustment arguments")

    # Box adjustment arguments
    box_adjustment.add_argument("--custom_box",
                                metavar="<file>",
                                type=str,
                                help="Custom search box based on a tab-separated file containing the minimum and \
                                     maximum cartesian values of each axis in angstrom.")
    box_adjustment.add_argument("--residues_box",
                                metavar="<file>",
                                type=str,
                                help="Custom search box around a list of residues (resnum_chainId) with a padding, \
                                     which are defined through a tab-separated file.")
    box_adjustment.add_argument("--padding",
                                metavar="<float>",
                                default=3.5,
                                type=positive_float,
                                help=u"Length (\u212B) added in each box direction. (default: %(default).1lf)")

    # Create argument group
    ligand_adjustment = parser.add_argument_group("Ligand adjustment arguments")

    # Ligand adjustment arguments
    ligand_adjustment.add_argument("--ligand",
                                   metavar="<.pdb>",
                                   type=str,
                                   help="Target-ligand trajectory to limit a search space within a radius \
                                        (ligand_cutoff) around it.")
    ligand_adjustment.add_argument("--ligand_cutoff",
                                   metavar="<float>",
                                   default=5.0,
                                   type=positive_float,
                                   help="Radius used to limit search space around the target-ligand. \
                                        (default: %(default).1lf)")

    return parser