"""
Plot your lpath results.
"""
import logging
import numpy
import lpath
import matplotlib
import matplotlib.pyplot as plt
from lpath.match import gen_dist_matrix, hcluster
from lpath.io import load_file
from os.path import exists

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

log = logging.getLogger(__name__)


def plt_config(args):
    if exists(args.mpl_styles):
        plt.style.use(args.mpl_styles)
    else:
        style_f = (files(lpath) / 'data/styles/default.mplstyle')
        plt.style.use(style_f)
        log.debug(f'DEBUG: Using default {style_f}')

    fig, ax = plt.subplots()

    return fig, ax


def plot(data, arguments):
    """
    This is an example function for plotting things.

    Parameters
    ----------
    input_array : numpy.ndarray or list
        An array generated from expanded_load.

    Returns
    -------
    """
    fig, ax = plt_config(arguments)

    return fig, ax


def organize_data(arguments):
    cluster_labels = load_file(arguments.cl_output, 1)
    load_data = load_file(arguments.output_pickle, 1)

    data = dict()

    npathways = len(load_data)
    log.info(f'There are {npathways} pathways.')
    path_idxs = numpy.arange(0, npathways)

    weights, durations, pathways = [], [], []
    for pathway in load_data:
        weights.append(pathway[0][-1])
        durations.append(int(len(pathway[:])))
        pathways.append(pathway)

    data['weights'] = numpy.asarray(weights)
    data['durations'] = numpy.asarray(durations)
    data['pathways'] = numpy.asarray(pathways)

    return data


def main(arguments):
    """
    Main function that executes the whole ``plot`` step.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """

    data = organize_data(arguments)

    plot(data, arguments)

    pass
