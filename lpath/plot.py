"""
Plot your lpath results.
"""
import logging
import numpy
import lpath
import matplotlib
import matplotlib.pyplot as plt
from lpath.match import determine_reassign, process_shorter_traj, gen_dist_matrix, hcluster
from lpath.io import load_file
from os.path import exists

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

log = logging.getLogger(__name__)


def process_plot_args(arguments):
    """
    Process plot arguments.

    Parameters
    ----------
    arguments : argparse.Namespace
        Arguments Namespace object as processed by argparser.

    Returns
    -------
    reassign : function
        The reassign function.

    """
    reassign = determine_reassign(arguments.plot_reassign_method)

    # In cases where you're not remaking the distance matrix but didn't provide a file path.
    if arguments.dmatrix_save is None:
        setattr(arguments, 'dmatrix_save', 'succ_traj/distmat.npy')

    return reassign


def plt_config(arguments):
    """
    Process matplotlib arguments and create fig/axis objects.

    Parameters
    ----------
    arguments : argparse.Namespace
        Arguments Namespace object as processed by argparser.

    Returns
    -------
    fig : matplotlib.pyplot.figure
        The matplotlib.pyplot.figure object for plotting.

    ax : matplotlib.axes
        The matplotlib.axes object for plotting.

    """
    if exists(arguments.mpl_styles):
        plt.style.use(arguments.mpl_styles)
    else:
        style_f = (files(lpath) / 'data/styles/default.mplstyle')
        plt.style.use(style_f)
        log.debug(f'DEBUG: Using default {style_f}')

    fig, ax = plt.subplots()

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

    return data, pathways


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


def main(arguments):
    """
    Main function that executes the whole ``plot`` step.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    # I/O and stuff...
    data, pathways = organize_data(arguments)
    reassign = process_plot_args(arguments)

    # Reassignment... (or not).
    dictionary = {}
    dictionary = reassign(data, pathways, dictionary, arguments.assign_name)  # system-specific reassignment of states

    if len(dictionary) < 3:
        log.warning(f'WARNING: Only {len(dictionary)} states defined, including the "unknown" state. \
                      This will likely produce bad clustering results and you should considering reassigning to more \
                      intermediate states using a modified ``--reassign-method``.')

    log.debug(f'Completed reassignment.')

    # Cleanup
    test_obj = process_shorter_traj(pathways, dictionary, arguments.exclude_short)  # Necessary if pathways are of variable length
    log.debug(f'Cleaned up trajectories.')

    if arguments.plot_remove_ends:
        test_obj = numpy.asarray([i[1:-1] for i in test_obj])

    dist_matrix, weights = gen_dist_matrix(test_obj, dictionary, file_name=arguments.dmatrix_save,
                                           out_dir=arguments.out_dir,
                                           remake=arguments.plot_dmatrix_remake,  # Calculate distance matrix
                                           metric=arguments.plot_longest_subsequence,  # Which metric to use
                                           condense=arguments.plot_condense,  # Whether to condense state strings
                                           n_jobs=arguments.plot_dmatrix_parallel)  # Number of parallel jobs


    fig, ax = plot(data, arguments)

    pass
