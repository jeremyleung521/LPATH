"""
Plot your lpath results.
"""
import logging
import numpy
import lpath
import matplotlib
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from lpath.io import load_file
from lpath.match import gen_dist_matrix, visualize, determine_rerun, ask_number_cluster, hcluster
from os.path import exists

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

log = logging.getLogger(__name__)


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

    return pathways


def determine_relabel(relabel_method):
    """
    Process argparse arguments to determine function to relabel trajectories.

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


class LPATHPlot:
    """
    A class consisting of a bunch of pre-made data for plotting.

    """

    def __init__(self, arguments):
        """
        Loads in and organizes data...

        """
        # Loading in things
        self.pathways = load_file(arguments.output_pickle)
        loaded_dmatrix = load_file(arguments.dmatrix_save) or None
        if isinstance(loaded_dmatrix, (numpy.ndarray, list)):
            # Turn into 1-D condensed distance matrix
            self.dist_matrix = squareform(loaded_dmatrix, checks=False) if loaded_dmatrix.ndim == 2 else loaded_dmatrix
        else:
            log.warning(f'No distance matrix inputted. Will proceed to rebuild.')
            self.dist_matrix = None
        self.cluster_labels = load_file(arguments.cl_output) if not arguments.regen_cl else None

        # Preparing empty array variables
        weights, durations, path_indices, iter_num = [], [], [], []

        # Report number of pathways.
        self.n_pathways = len(self.pathways)
        log.info(f'There are {self.n_pathways} pathways.')

        for pathway in loaded_data:
            weights.append(pathway[0][-1])
            durations.append(len(pathway[:]))
            iter_num += [frame[0] for frame in pathway]
        iter_num = numpy.asarray(iter_num)

        for cluster in set(self.cluster_labels):
            path_indices.append(numpy.where(self.cluster_labels == cluster)[0])

        # Precalculating stuff for later...
        self.arguments = arguments
        self.out_path = arguments.out_path
        self.weights = numpy.asarray(weights)
        self.durations = numpy.asarray(durations)
        self.path_indices = path_indices

        self.min_iter = min(iter_num[numpy.nonzero(iter_num)])
        self.max_iter = max(iter_num)
        self.interval = 10 if self.max_iter - self.min_iter > 50 else 2
        self.dendrogram_threshold = arguments.dendrogram_threshold
        self.new_pathways = self.pathways
        self.mpl_colors = None if arguments.mpl_colors is ['None'] else arguments.mpl_colors

        # Set up dummy variables for the plots!
        self.fig = None
        self.ax = None
        self.show_fig = arguments.dendrogram_show

    def plt_config(self):
        """
        Process matplotlib arguments and append fig/axis objects to class.

        Returns
        -------
        self.fig : matplotlib.Figure
            The matplotlib.Figure object for plotting.

        self.ax : matplotlib.axes
            The matplotlib.Axes object for plotting.

        """
        if exists(self.arguments.mpl_styles):
            plt.style.use(self.arguments.mpl_styles)
        else:
            style_f = (files(lpath) / 'data/styles/default.mplstyle')
            plt.style.use(style_f)
            log.debug(f'DEBUG: Using default {style_f}')

        self.fig, self.ax = plt.subplots()

    def plot(self):
        """
        This is an example method for plotting things.

        """
        self.plt_config()

    def plot_branch_colors(self):
        if self.fig is None or self.ax is None:
            self.plt_config()

        z = sch.linkage(self.dist_matrix, method='ward')

        sch.set_link_color_palette(self.mpl_colors)
        # cluster_colors_array = [self.mpl_colors[label] for label in self.cluster_labels]
        # link_cols = {}

        # for i, i_pair in enumerate(z[:, :2].astype(int)):
        #    c1, c2 = (link_cols[x] if x > len(z) else cluster_colors_array[x] for x in i_pair)
        #    link_cols[i + 1 + len(z)] = c1 if c1 == c2 else self.mpl_colors[-1]

        try:
            # Temporarily override the default line width:
            with plt.rc_context({'lines.linewidth': 3}):
                sch.dendrogram(z, no_labels=True, color_threshold=self.dendrogram_threshold,
                               above_threshold_color=self.mpl_colors[-1], ax=self.ax)
        except RecursionError as e:
            # Catch cases where are too many branches in the dendrogram for default recursion to work.
            import sys

            sys.setrecursionlimit(100000)
            log.warning(e)
            log.warning(f'WARNING: Dendrogram too complex to plot with default settings. Upping the recursion limit.')
            # Temporarily override the default line width:
            with plt.rc_context({'lines.linewidth': 3}):
                sch.dendrogram(z, no_labels=True, color_threshold=self.dendrogram_threshold,
                               above_threshold_color=self.mpl_colors[0], ax=self.ax)

        self.ax.axhline(y=self.dendrogram_threshold, c='k', linestyle='--', linewidth=2.5)
        self.ax.set_ylabel("distance")
        self.ax.set_xlabel("pathways")
        self.fig.set_layout_engine(layout='tight')
        self.fig.savefig(f'{self.out_path}/dendrogram_custom_color.pdf')

        if self.show_fig:
            plt.show()


def plot_custom(fig, ax):
    abs_pcoord = numpy.abs(numpy.diff(pcoord))
    mask = numpy.hstack([abs_pcoord > 250, [False]])
    output[jdx] = numpy.ma.MaskedArray(pcoord1, mask)

    plt.plot(masked_pcoord1, masked_pcoord2, color=colors[cidx], alpha=0.25)


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
        log.debug('Setting distance matrix output to default.')

    if arguments.cl_output is None:
        setattr(arguments, 'cl_output', 'succ_traj/cluster_labels.npy')
        log.debug('Setting cluster label output to default.')

    if arguments.output_pickle is None:
        setattr(arguments, 'output_pickle', 'succ_traj/pathways.pickle')
        log.debug('Setting match pickle output/plot pickle input to default.')

    if arguments.dendrogram_threshold is None:
        setattr(arguments, 'dendrogram_threshold', 0.5)
        log.debug('Setting dendrogram threshold output to default 0.5.')

    return relabel


def plot_hist_iter_num(data, fig, ax, mpl_colors, out_path):
    """

    Parameters
    ----------
    data : LPATHPlot
        A PlotData class containing a bunch of information for plotting.

    fig : matplotlib.Figure
        The matplotlib.Figure object for plotting.

    ax : matplotlib.axes
        The matplotlib.Axes object for plotting.

    mpl_colors : list of colors
        A list of colors for plotting. The last color is typically reserved for special graying out.

    out_path : str, default: 'graphs'
        Path to output the files.

    """
    for cluster in range(len(data.path_indices)):
        iteration_target = []
        path_idxs_c = data.pathways[data.path_indices[cluster]]
        for idx, pathway in enumerate(data.pathways):
            if idx in path_idxs_c:
                iteration_target.append(pathway[0])
        ax.hist(iteration_target, bins=numpy.arange(data.min_iter - 1, data.max_iter, data.interval),
                weights=data.weights[data.path_indices], color=mpl_colors[cluster], alpha=0.7, ax=ax)

    ax.set_xlim(data.min_iter - 1, data.max_iter)
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
    relabel = process_plot_args(arguments)
    data = LPATHPlot(arguments)

    # Reassignment... (or not).
    data.new_pathways, data.dictionary = relabel(data.pathways)  # system-specific relabeling of things

    if arguments.regen_cl:
        # data.dist_matrix, data.weights = gen_dist_matrix(data.new_pathways, data.dictionary,
        #                                                  file_name=arguments.dmatrix_save,
        #                                                  out_dir=arguments.out_dir,
        #                                                  remake=True,  # Calculate distance matrix
        #                                                  metric=arguments.plot_longest_subsequence,  # Which metric to use
        #                                                  condense=arguments.plot_condense,
        #                                                  # Whether to condense consecutive state strings
        #                                                  n_jobs=arguments.dmatrix_parallel,
        #                                                  # Number of jobs for pairwise_distance
        #                                                  )

        # Attempt to visualize dendrogram with new the distance matrix.
        data.ax = visualize(data.dist_matrix, threshold=arguments.dendrogram_threshold, out_dir=arguments.out_dir,
                            show_fig=arguments.dendrogram_show, mpl_colors=arguments.mpl_colors, ax=data.ax)
        determine_rerun(data.dist_matrix, mpl_colors=arguments.mpl_colors, ax=data.ax)
        ncluster = ask_number_cluster()
        cluster_labels = hcluster(data.dist_matrix, ncluster)
        data.cluster_labels = cluster_labels

    # Cleanup
    # test_obj = process_shorter_traj(pathways, dictionary, arguments.plot_exclude_short)
    # Necessary if pathways are of variable length
    # log.debug(f'Cleaned up trajectories.')

    data.plot_branch_colors()
