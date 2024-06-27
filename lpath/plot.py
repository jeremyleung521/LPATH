"""
Plot your lpath results.
"""
import numpy
import lpath
import matplotlib
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import networkx
from scipy.spatial.distance import squareform
from os.path import exists
from ast import literal_eval

from lpath.extloader import *
from lpath.io import load_file
from lpath.match import calc_linkage, visualize, determine_rerun, ask_number_clusters, hcluster

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

from ._logger import Logger

log = Logger().get_logger(__name__)


def relabel_identity(data):
    """
    Use ``lpath.match`` states as is. Does not attempt to relabel anything.

    Parameters
    ----------
    data : LPATHPlot class
        An LPATHPlot class object.

    Returns
    -------
    pathways : numpy.ndarray
        A (not-modified) array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    cluster_labels : numpy.ndarray
        A (not-modified) array of cluster labels. Passed here just in case you want to do something fancy.

    """
    return data.pathways, data.cluster_labels


def relabel_custom(data):
    """
    Relabel pathways (pcoord or states) from ``lpath match`` frames into different values. This is highly
    specific to the system. If ``lpath match``'s definition is sufficient,
    you can proceed with what's made in the previous step using ``relabel_identity``.

    In this example, we're modifying it so the phi/psi angles (columns 3 and 4) are in (-180,180] instead.

    Parameters
    ----------
    data : LPATHPlot class
        An LPATHPlot class object.

    Returns
    -------
    pathways : numpy.ndarray
        A modified array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    cluster_labels : numpy.ndarray
        A modified array of cluster labels. Passed here just in case you want to do something fancy.

    """
    # Looping through pathways and modifying it
    for pathway in data.pathways:
        # For each pathway
        for pcoord in [pathway[:, 3], pathway[:, 4]]:
            # Changing forth and fifth column, which are the phi/psi values from pcoord.
            for pidx, val in enumerate(pcoord):
                # Looping through each frame
                if val > 180:
                    # If it's > 180, subtract 360
                    pcoord[pidx] -= 360

    return data.pathways, data.cluster_labels


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
        try:
            loaded_dmatrix = load_file(arguments.dmatrix_save)
        except (ValueError, FileNotFoundError):
            loaded_dmatrix = None
        if isinstance(loaded_dmatrix, (numpy.ndarray, list)):
            # Turn into 1-D condensed distance matrix
            self.dist_matrix = squareform(loaded_dmatrix, checks=False) if loaded_dmatrix.ndim == 2 else loaded_dmatrix
            self.linkage = calc_linkage(self.dist_matrix)
        else:
            raise ValueError(f'No distance matrix provided. Please generate one with ``lpath match [...]`` and '
                             'pass it to ``lpath plot`` with the ``--plot-input`` flag.')
        self.cluster_labels = load_file(arguments.cl_output)  # if not arguments.regen_cl else None

        # Preparing empty array variables
        weights, durations, path_indices, iter_num, target_iter, states = [], [], [], [], [], []

        # Report number of pathways.
        self.n_pathways = len(self.pathways)
        log.info(f'There are {self.n_pathways} pathways.')

        for pathway in self.pathways:
            non_zero = pathway[numpy.nonzero(pathway[:, 0])]  # Removing padding frames...
            weights.append(non_zero[-1, -1])
            durations.append(len(non_zero) / arguments.stride)
            target_iter.append(non_zero[-1][0])
            iter_num += [frame[0] for frame in non_zero]
            states.append(pathway[:, 2])

        for cluster in set(self.cluster_labels):
            path_indices.append(numpy.where(self.cluster_labels == cluster)[0])

        # Precalculating stuff for later...
        self.arguments = arguments
        self.out_path = arguments.out_path
        self.weights = numpy.asarray(weights)
        self.durations = numpy.asarray(durations)
        self.states = numpy.asarray(states, dtype=object)
        self.path_indices = path_indices
        self.num_iter = numpy.asarray(iter_num, dtype=object)
        self.target_iter = numpy.asarray(target_iter)
        self.num_clusters = len(path_indices)
        self.stride = arguments.stride

        # weighted_counts = [numpy.sum(self.weights[self.states == i]) for i in range(self.num_clusters)]
        # self.weighted_counts = numpy.array(weighted_counts, dtype=float)

        self.min_iter = min(self.num_iter)
        self.max_iter = max(self.num_iter)
        self.interval = 25 if self.max_iter - self.min_iter > 50 else 5
        self.bins = int((self.max_iter - self.min_iter) // self.interval)
        self.dendrogram_threshold = arguments.dendrogram_threshold
        self.new_pathways = self.pathways
        self.mpl_colors = None if arguments.mpl_colors is ['None'] else arguments.mpl_colors

        try:
            self.mpl_args = literal_eval(arguments.matplotlib_args)
        except SyntaxError:
            self.mpl_args = dict()

        # Set up dummy variables for the plots!
        self.fig = None
        self.ax = None
        self.show_fig = arguments.dendrogram_show
        self.ax_idx = 0

        if not self.show_fig:
            matplotlib.use('agg')

    def plt_config(self, ax_idx=None, separate=None):
        """
        Process matplotlib arguments and append fig/axis objects to class.

        ax_idx = list or int
            Index or list of indic

        """
        if exists(self.arguments.mpl_styles):
            plt.style.use(self.arguments.mpl_styles)
        else:
            style_f = (files(lpath) / 'data/styles/default.mplstyle')
            plt.style.use(style_f)
            log.info(f'Using default {style_f}')

        if separate:
            try:
                mpl_keys = self.arguments.matplotlib_args.keys()
            except AttributeError:
                mpl_keys = {}.keys()

            if 'nrows' in mpl_keys and 'ncols' in mpl_keys:
                log.info(f'Defaulting to user-provided plot layout.')
                n_rows = self.arguments.matplotlib_args['nrows']
                n_cols = self.arguments.matplotlib_args['ncols']
            else:
                log.info(f'Defaulting to default LPATH plot layout')
                if self.num_clusters > 4:
                    n_cols = 3
                else:
                    n_cols = 2

                n_rows = -(-self.num_clusters // n_cols)

            # Doing ceiling division to determine how many rows available.
            self.fig, self.ax = plt.subplots(nrows=n_rows, ncols=n_cols)
        else:
            self.fig, self.ax = plt.subplots(**self.mpl_args)

        if isinstance(self.ax, plt.Axes):
            # If only a single axes, put into list.
            self.ax = numpy.asarray([self.ax], dtype=object)

        self.ax = self.ax.flatten()

    def determine_plot_axes(self, ax_idx=None, separate=False):
        """
        Determine which axes to plot and return a list of axes to plot.

        Parameter
        ---------
        ax_idx : list of int or None, default: None
            Which axes index to plot the graph.

        separate : bool, default: False
            Whether to plot each cluster in separate subplots or not.

        """
        if separate is False:
            if ax_idx:
                plot_axes = [self.ax[ax_idx[0]]]

            else:
                # Plotting all into first axes.
                plot_axes = [self.ax[self.ax_idx]]

        else:
            # Plot all into subsequent axes
            if ax_idx is not None:
                if len(ax_idx) >= len(self.path_indices):
                    plot_axes = self.ax[ax_idx]
                    self.ax_idx += len(ax_idx) - 1
                else:
                    log.warning(f'Not enough axes to plot each one individually. Plotting all into first idx provided.')
                    plot_axes = [self.ax[ax_idx[self.ax_idx]]]
                    self.ax_idx += 1

            else:
                log.warning(f'Plot axes index not given. Plotting all into first N axes.')
                plot_axes = self.ax[list(range(self.num_clusters))]
                self.ax_idx += 1

        return plot_axes

    def plot(self):
        """
        This is an example method for plotting things. You can set up subplots with
        fig and ax first with plt_config().

        """
        self.plt_config()

        self.ax.plot(self.pathways[self.path_indices[0], 3], self.pathways[self.path_indices[0], 3])
        self.fig.savefig('output.pdf', dpi=300)

    def plotdendro_branch_colors(self, ax_idx=None):
        """
        Plot dendrogram branches with customized colors defined by ``--mpl_colors``.

        Parameters
        ----------
        ax_idx : list of int or None, default: None
            Which axes index to plot the graph.

        """
        self.plt_config(ax_idx)

        # Code from stackoverflow to recolor each branch and node with specific colors.
        # Setting link color palette is much simpler.
        # cluster_colors_array = [self.mpl_colors[label] for label in self.cluster_labels]
        # link_cols = {}
        # len_link = len(self.linkage)
        # for i, i_pair in enumerate(self.linkage[:, :2].astype(int)):
        #    c1, c2 = (link_cols[x] if x > len_link else cluster_colors_array[x] for x in i_pair)
        #    link_cols[i + 1 + len(self.linkage)] = c1 if c1 == c2 else self.mpl_colors[-1]
        sch.set_link_color_palette(self.mpl_colors[:-1])

        plot_axes = self.determine_plot_axes(ax_idx, False)

        try:
            # Temporarily override the default line width:
            with plt.rc_context({'lines.linewidth': 2}):
                sch.dendrogram(self.linkage, no_labels=True, color_threshold=self.dendrogram_threshold,
                               above_threshold_color=self.mpl_colors[-1], ax=plot_axes[0])
        except RecursionError as e:
            # Catch cases where are too many branches in the dendrogram for default recursion to work.
            import sys

            sys.setrecursionlimit(100000)
            log.warning(e)
            log.warning(f'WARNING: Dendrogram too complex to plot with default settings. Upping the recursion limit.')
            # Temporarily override the default line width:
            with plt.rc_context({'lines.linewidth': 2}):
                sch.dendrogram(self.linkage, no_labels=True, color_threshold=self.dendrogram_threshold,
                               above_threshold_color=self.mpl_colors[-1], ax=plot_axes[0])

        plot_axes[0].axhline(y=self.dendrogram_threshold, c='k', linestyle='--', linewidth=1.5)
        plot_axes[0].set_ylabel("distance")
        plot_axes[0].set_xlabel("pathways")
        self.fig.set_layout_engine(layout='tight')
        self.fig.savefig(f'{self.out_path}/dendrogram.pdf')

        if self.show_fig:
            plt.show()

    def plothist_weight_cluster(self, ax_idx=None):
        """
        Plot histogram of cluster number vs. weight.

        Parameter
        ---------
        ax_idx : list of int or None, default: None
            Which axes index to plot the graph.

        """
        self.plt_config(ax_idx)
        plot_axes = self.determine_plot_axes(ax_idx)

        for axes in plot_axes:
            x_list = []
            y_list = []
            for c_idx, cluster_indices in enumerate(self.path_indices):
                x_list.append(c_idx + 1)
                y_list.append(numpy.sum(self.weights[cluster_indices]))

            axes.bar(x_list, y_list, width=0.75, color=self.mpl_colors[:-1])
            axes.set_xlim(0, max(x_list) + 1)
            # axes.set_ylim(0, max(y_list))
            axes.set_xticks(ticks=x_list, labels=[f'{x}' for x in x_list])
            axes.set_xlabel(r'pathway classes')
            axes.set_ylabel('probability')

        self.fig.set_layout_engine(layout='tight')
        self.fig.savefig(f"{self.out_path}/weight_histogram.pdf", dpi=300)
        log.info(f'Outputted graph in {self.out_path}/weight_histogram.pdf.')

        if self.show_fig:
            plt.show()

    def plothist_event_duration(self, ax_idx=None, separate=False):
        """
        Plot histogram of vent duration time vs. iteration number history number
        with customized colors defined by ``--mpl_colors``.

        Parameter
        ---------
        ax_idx : list of int or None, default: None
            Which axes index to plot the graph.

        separate : bool, default: False
            Whether to plot each cluster in separate subplots or not.

        """
        self.plt_config(ax_idx, separate=True)
        plot_axes = self.determine_plot_axes(ax_idx, separate)

        for n_axes, axes in enumerate(plot_axes):
            if separate:
                # Plot each pathway class in a separate axes.
                axes.hist(self.durations[self.path_indices[n_axes]], bins=self.bins,
                          range=(self.min_iter, self.max_iter),
                          weights=self.weights[self.path_indices[n_axes]], density=True,
                          color=self.mpl_colors[n_axes % (len(self.mpl_colors) - 1)])
                axes.set_title(f'class {n_axes}')
            else:
                # Plot all pathway classes into same axes
                axes.hist(self.durations[c_indices], bins=self.bins,
                          weights=self.weights[c_indices], density=True,
                          color=self.mpl_colors[cluster_id % (len(self.mpl_colors) - 1)])
            axes.set_xlim(0, max(self.durations))
            # axes.set_ylim(0, max(self.weights))
            axes.set_xlabel('duration')
            axes.set_ylabel('probability')

        self.fig.set_layout_engine(layout='tight')
        self.fig.savefig(f"{self.out_path}/durations.pdf", dpi=300)
        log.info(f'Outputted graph in {self.out_path}/durations.pdf.')

        if self.show_fig:
            plt.show()

    def plothist_target_iter(self, ax_idx=None, separate=False):
        """
        Plot network of transitions between different states.

        Parameter
        ---------
        ax_idx : list of int or None, default: None
            Which axes index to plot the graph.

        separate : bool, default: False
            Whether to plot each cluster in separate subplots or not.

        """
        self.plt_config(ax_idx, separate=True)
        plot_axes = self.determine_plot_axes(ax_idx, separate)

        for n_axes, axes in enumerate(plot_axes):
            iteration_target = self.target_iter[self.path_indices[n_axes]]
            axes.hist(iteration_target, bins=numpy.arange(self.min_iter - 1, self.max_iter, self.interval),
                      weights=self.weights[self.path_indices[n_axes]],
                      color=self.mpl_colors[n_axes % (len(self.mpl_colors) - 1)],
                      label=f'class {n_axes}')

            axes.set_xlim(self.min_iter - 1, self.max_iter)
            axes.set_xlabel("iteration")
            axes.set_ylabel("probability")
            axes.set_yscale("log")

        self.fig.set_layout_engine(layout='tight')
        self.fig.savefig(f"{self.out_path}/target_iteration.pdf")
        log.info(f'Outputted graph in {self.out_path}/target_iteration.pdf.')

        if self.show_fig:
            plt.show()

    def plot_network(self, ax_idx=None, separate=False):
        """
        Plot network of target iteration vs. iteration number history number
        with customized colors defined by ``--mpl_colors``.

        Parameter
        ---------
        ax_idx : list of int or None, default: None
            Which axes index to plot the graph.

        separate : bool, default: False
            Whether to plot each cluster in separate subplots or not.

        """
        self.plt_config(ax_idx, separate=True)
        plot_axes = self.determine_plot_axes(ax_idx, separate)

        for n_axes, axes in enumerate(plot_axes):
            G = networkx.DiGraph()
            colors = [
                matplotlib.colors.cnames[default_dendrogram_colors[cidx]] if mapping[i] not in ["A", "B"] else "#D3D3D3"
                for i in range(n_total_states)]

            iteration_target = self.target_iter[self.path_indices[n_axes]]
            axes.hist(iteration_target, bins=numpy.arange(self.min_iter - 1, self.max_iter, self.interval),
                      weights=self.weights[self.path_indices[n_axes]],
                      color=self.mpl_colors[n_axes % (len(self.mpl_colors) - 1)],
                      label=f'class {n_axes}')

            axes.set_xlim(self.min_iter - 1, self.max_iter)
            axes.set_xlabel("iteration")
            axes.set_ylabel("probability")
            axes.set_yscale("log")

        self.fig.set_layout_engine(layout='tight')
        self.fig.savefig(f"{self.out_path}/network.pdf")
        log.info(f'Outputted graph in {self.out_path}/network.pdf.')

        if self.show_fig:
            plt.show()


def plot_custom():
    """
    Example custom function for custom plot script for plotting with the LPATHPlot class.
    In here, we plot each cluster in a separate subplot, and pathway onto a phi

    """
    import argparse

    args = argparse.Namespace(
        output_pickle='succ_traj/pathways.pickle',  # Output pickle object from lpath.match
        cl_output='succ_traj/cluster_labels.npy',  # Cluster label outputs from lpath.match
        dmatrix_save='succ_traj/distmat.npy',  # Distance matrix output from lpath.match
        matplotlib_args={'nrows': 1, 'ncols': 2},  # Arguments for subplots
    )

    data = LPATHPlot(args)

    data.plt_config()

    # Example for
    n_states = int(max([seg[2] for traj in data.pathways for seg in traj])) + 1  # Plus 1 for 0-indexing

    # Looping through
    for cluster in range(n_states):
        # Looping through each cluster
        for p_idx in data.path_indices[cluster]:
            # Looping through each pathway
            output = []
            for dimension in [3, 4]:
                # Plotting both Phi / Psi, which is columns 3, 4 in the pickle object from lpath.match
                pcoord = data.pathways[p_idx][:, dimension]
                abs_pcoord = numpy.abs(numpy.diff(pcoord[dimension]))
                mask = numpy.hstack([abs_pcoord > 250, [False]])
                output.append(numpy.ma.MaskedArray(pcoord, mask))
            data.ax[cluster].plot(output[0], output[1], color=data.mpl_colors[cluster % len(data.mpl_colors)],
                                  alpha=0.25)

    data.fig.savefig('custom_graph.pdf', dpi=300)


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
        The relabeling function.

    """
    relabel = determine_relabel(arguments.relabel_method)

    # In cases where you're not going through ``match`` first, some arguments might be empty.
    if arguments.dmatrix_save is None:
        setattr(arguments, 'dmatrix_save', 'succ_traj/distmat.npy')
        log.warning(f'Setting distance matrix output to default {arguments.dmatrix_save}. Make sure you\'re sure '
                    'of this, or remake the distance matrix with ``lpath match``.')

    if arguments.cl_output is None:
        setattr(arguments, 'cl_output', 'succ_traj/cluster_labels.npy')
        log.info(f'Setting cluster label output to default {arguments.cl_output}.')

    if arguments.output_pickle is None:
        setattr(arguments, 'output_pickle', 'succ_traj/pathways.pickle')
        log.info(f'Setting match pickle output/plot pickle input to default {arguments.output_pickle}.')

    if arguments.dendrogram_threshold is None:
        setattr(arguments, 'dendrogram_threshold', 0.5)
        log.info(f'Setting dendrogram threshold output to default {arguments.dendrogram_threshold}.')

    return relabel


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
    data.new_pathways, data.new_cluster_labels = relabel(data)  # system-specific relabeling of things

    # if arguments.regen_cl:
    #     # Attempt to re-visualize the dendrogram.
    #     data.ax = visualize(data.linkage, threshold=arguments.dendrogram_threshold, out_path=arguments.out_path,
    #                         show_fig=arguments.dendrogram_show, mpl_colors=arguments.mpl_colors, ax=data.ax)
    #     determine_rerun(data.linkage, out_path=data.out_path, mpl_colors=arguments.mpl_colors, ax=data.ax)
    #     data.num_clusters = ask_number_clusters(arguments.num_clusters)
    #     cluster_labels = hcluster(data.linkage, data.num_clusters)
    #     data.cluster_labels = cluster_labels
    #     numpy.save(f'{data.out_path}/new_cluster_labels.npy', data.cluster_labels)

    # The following plots are done in a single plot
    data.plotdendro_branch_colors()
    data.plothist_weight_cluster()

    # The following plots are done per-cluster
    data.plothist_target_iter(separate=True)
    data.plothist_event_duration(separate=True)

