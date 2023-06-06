"""
Plot your lpath results.
"""
import logging
import numpy
import lpath
import matplotlib
import matplotlib.pyplot as plt
from lpath.io import load_file
from lpath.match import visualize, determine_rerun, ask_number_cluster, hcluster
from os.path import exists

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

log = logging.getLogger(__name__)

default_dendrogram_colors = ['tomato', 'dodgerblue', 'red', 'purple', 'grey']


class plotData:
    """
    A class consisting of a bunch of pre-made data for plotting.

    """
    def __init__(self, arguments):
        """
        Loads in and organizes data...

        """
        # Loading in things
        cluster_labels = load_file(arguments.cl_output)
        loaded_data = load_file(arguments.output_pickle)
        loaded_dmatrix = load_file(arguments.dmatrix_save)

        # Preparing variables
        data = dict()
        weights, durations, path_indices, iter_num = [], [], [], []

        # Report pathways.
        n_pathways = len(loaded_data)
        log.info(f'There are {n_pathways} pathways.')

        for pathway in loaded_data:
            weights.append(pathway[0][-1])
            durations.append(len(pathway[:]))
            iter_num += [frame[0] for frame in pathway]
        iter_num = numpy.asarray(iter_num)

        for cluster in set(cluster_labels):
            path_indices.append(numpy.where(cluster_labels == cluster)[0])

        # Precalculating stuff for later...
        self.weights = numpy.asarray(weights)
        self.durations = numpy.asarray(durations)
        self.pathways = loaded_data
        self.path_indices = path_indices
        self.dist_matrix = loaded_dmatrix
        self.cluster_labels = cluster_labels
        self.n_pathways = len(loaded_data)
        self.min_max_iter = [min(iter_num[numpy.nonzero(iter_num)]), max(iter_num)]



def relabel_identity(pathways, cluster_labels):
    """
    Use ``lpath.match`` states as is. Does not attempt to relabel anything.

    Parameters
    ----------
    pathways : numpy.ndarray
        An array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    Returns
    -------
    pathways : numpy.ndarray
        A modified array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    """
    return pathways


def relabel_custom(pathways, cluster_labels):
    """
    Relabel pathways (pcoord or states) from ``lpath.match`` frames into different values. This is highly
    specific to the system. If ``lpath.match``'s definition is sufficient,
    you can proceed with what's made in the previous step using ``relabel_identity``.

    In this example, we're modifying it so the phi/psi angles (columns 3 and 4) are in (-180,180] instead.

    Parameters
    ----------
    pathways : numpy.ndarray
        An array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    Returns
    -------
    pathways : numpy.ndarray
        A modified array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    """
    # Example for relabeling pcoord and generating the mapping dictionary.
    n_states = int(max([seg[2] for traj in pathways for seg in traj])) + 1

    # Looping through
    for cidx, cluster in range(n_states):
        path_idxs_c = path_idxs[cluster_labels == cluster]
        for idx, pathway in enumerate(pathways):
            if idx in path_idxs_c:
                pathway = numpy.array(pathway)
                pcoord1 = pathway[:, 3]
                pcoord2 = pathway[:, 4]

                for jdx, pcoord in enumerate([pcoord1, pcoord2]):
                    for pidx, val in enumerate(pcoord):
                        if val > 180:
                            pcoord[pidx] -= 360

    dictionary = {}
    for idx in range(n_states):
        dictionary[idx] = str(idx)

    dictionary[n_states] = '!'  # Unknown state

    return pathways


def plot_custom(fig, ax):
    abs_pcoord = numpy.abs(numpy.diff(pcoord))
    mask = numpy.hstack([abs_pcoord > 250, [False]])
    output[jdx] = numpy.ma.MaskedArray(pcoord1, mask)

    plt.plot(masked_pcoord1, masked_pcoord2, color=colors[cidx], alpha=0.25)


def determine_relabel(relabel_method):
    """
    Argument processing to determine function to relabel trajectories.

    Parameters
    ----------
    relabel_method : str , default: 'relabel_identity'
        String from argument.plot_relabel_identity, straight from argparser.

    Returns
    -------
    relabel : function
        The relabelling function.

    """
    preset_relabel = {
        'relabel_identity': relabel_identity,
        'relabel_custom': relabel_custom,
    }

    if relabel_method in preset_relabel.keys():
        relabel = preset_relabel[relabel_method]
    else:
        import sys
        import os
        sys.path.append(os.getcwd())

        relabel = get_object(relabel_method)
        log.info(f'INFO: Replaced relabel() with {relabel_method}')

    return relabel


def process_plot_args(arguments):
    """
    Process plot arguments.

    Parameters
    ----------
    arguments : argparse.Namespace
        Arguments Namespace object as processed by argparser.

    Returns
    -------
    relabel : function
        The relabel function.

    """
    relabel = determine_relabel(arguments.plot_relabel_method)

    # In cases where you're not going through `match` first, some arguments might be empty.
    if arguments.dmatrix_save is None:
        setattr(arguments, 'dmatrix_save', 'succ_traj/distmat.npy')

    if arguments.cl_output is None:
        setattr(arguments, 'cl_output', 'succ_traj/cluster_labels.npy')

    if arguments.output_pickle is None:
        setattr(arguments, 'output_pickle', 'succ_traj/pathways.pickle')

    return relabel


def plt_config(arguments):
    """
    Process matplotlib arguments and create fig/axis objects.

    Parameters
    ----------
    arguments : argparse.Namespace
        Arguments Namespace object as processed by argparser.

    Returns
    -------
    fig : matplotlib.Figure
        The matplotlib.Figure object for plotting.

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
    """
    Catch all function for organizing the data loaded in for plotting.

    Parameters
    ----------
    arguments : argparse.Namespace
        An argparse.Namespace object straight from argparser.

    """

    data = plotData(load_file(arguments.cl_output), load_file(arguments.output_pickle), load_file(arguments.dmatrix_save))

    # Loading in things
    cluster_labels = load_file(arguments.cl_output)
    loaded_data = load_file(arguments.output_pickle)
    loaded_dmatrix = load_file(arguments.dmatrix_save)

    # Preparing variables
    data = dict()
    weights, durations, path_indices, iter_num = [], [], [], []

    # Report pathways.
    n_pathways = len(loaded_data)
    log.info(f'There are {n_pathways} pathways.')

    for pathway in loaded_data:
        weights.append(pathway[0][-1])
        durations.append(len(pathway[:]))
        iter_num += [frame[0] for frame in pathway]
    iter_num = numpy.asarray(iter_num)

    for cluster in set(cluster_labels):
        path_indices.append(numpy.where(cluster_labels == cluster)[0])

    # Precalculating stuff for later...
    data['weights'] = numpy.asarray(weights)
    data['durations'] = numpy.asarray(durations)
    data['pathways'] = loaded_data
    data['path_indices'] = path_indices
    data['dist_matrix'] = loaded_dmatrix
    data['cluster_labels'] = cluster_labels
    data['n_pathways'] = len(loaded_data)
    data['min_max_iter'] = [min(iter_num[numpy.nonzero(iter_num)]), max(iter_num)]

    return data


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


def plot_branch_colors(fig, ax, labels, mpl_colors=default_dendrogram_colors):
    z = sch.linkage(distmat_condensed, method="ward")
    cluster_colors_array = [mpl_colors[label] for label in labels]
    link_cols = {}
    for i, i12 in enumerate(z[:, :2].astype(int)):
        c1, c2 = (link_cols[x] if x > len(z) else cluster_colors_array[x] for x in i12)
        link_cols[i + 1 + len(z)] = c1 if c1 == c2 else mpl_colors[-1]

    try:
        # Temporarily override the default line width:
        with plt.rc_context({'lines.linewidth': 3}):
            sch.dendrogram(z, no_labels=True, color_threshold=threshold, link_color_func=lambda x: mpl_colors[x],
                           above_threshold_color=mpl_colors[-1], ax=ax)
    except RecursionError as e:
        # Catch cases where are too many branches in the dendrogram for default recursion to work.
        import sys

        sys.setrecursionlimit(100000)
        log.warning(e)
        log.warning(f'WARNING: Dendrogram too complex to plot with default settings. Upping the recursion limit.')
        with plt.rc_context({'lines.linewidth': 3}):
            sch.dendrogram(z, no_labels=True, color_threshold=threshold, link_color_func=lambda x: mpl_colors[x],
                           above_threshold_color=mpl_colors[-1], ax=ax)

    plt.axhline(y=threshold, c="k", linestyle="--", linewidth=2.5)
    plt.ylabel("distance")
    plt.xlabel("pathways")
    plt.tight_layout()
    plt.savefig("dendrogram_edit.pdf")


def plot_hist_iter_num(pathways, fig, ax, out_path):
    for cidx, cluster in enumerate([1, 2]):
        iteration_target = []
        path_idxs_c = path_idxs[cluster_labels == cluster]
        for idx, pathway in enumerate(pathways):
            if idx in path_idxs_c:
                pathway = numpy.array(pathway[0])
                iteration_target.append(pathway[0])
        ax.hist(iteration_target, bins=numpy.arange(0, 600, 10), weights=weights[cluster_labels == cluster],
                color=colors[cidx], alpha=0.7, ax=ax)

    ax.set_xlim(0, 600)
    ax.set_xlabel("we iteration of arrival")
    ax.set_ylabel("probability")
    ax.set_yscale("log")

    fig.savefig(f"{out_path}/iteration.pdf")
    log.info(f'Outputted Graph in {out_path}/iteration.pdf.')


def main(arguments):
    """
    Main function that executes the whole ``plot`` step.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    # I/O and stuff...
    data = organize_data(arguments)
    relabel = process_plot_args(arguments)

    # Reassignment... (or not).
    data['new_pathways'] = relabel(data['pathways'])  # system-specific relabeling of things

    if arguments.regen_cl:
        visualize(dist_matrix, threshold=arguments.dendrogram_threshold, out_dir=arguments.out_dir,
                  show=arguments.dendrogram_show, mpl_colors=arguments.mpl_colors)  # Visualize
        determine_rerun(dist_matrix, mpl_colors=arguments.mpl_colors)
        ncluster = ask_number_cluster()
        cluster_labels = hcluster(dist_matrix, ncluster)

    # if len(dictionary) < 3:
    #     log.warning(f'WARNING: Only {len(dictionary)} states defined, including the "unknown" state. \
    #                   This will likely produce bad clustering results and you should considering relabeling to more \
    #                   intermediate states using a modified ``--plot-relabel-method``.')

    log.debug(f'Completed relabeling.')

    # Cleanup
    # test_obj = process_shorter_traj(pathways, dictionary, arguments.plot_exclude_short)  # Necessary if pathways are of variable length
    # log.debug(f'Cleaned up trajectories.')

    fig, ax = plot(new_pathways)

    plot_branch_colors(fig, ax, cluster_labels, arguments.mpl_colors)

    pass
