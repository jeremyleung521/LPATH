"""
Plot your lpath results.
"""
import logging
import numpy
import lpath
import matplotlib
import matplotlib.pyplot as plt
from lpath.io import load_file
try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files
log = logging.getLogger(__name__)


def plot(input_array):
    """
    This is an example function for plotting things.

    Parameters
    ----------
    input_array : numpy.ndarray or list
        An array generated from expanded_load.

    Returns
    -------
    """
    inp_file = (files(lpath) / 'data/styles/default.mplstyle')
    print(inp_file)


def main(arguments):
    """
    Main function that executes the whole ``plot`` step.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    cluster_labels = load_file(arguments.cl_output, 1)
    pass


