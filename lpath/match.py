"""
Pattern match your extracted trajectories and cluster pathways classes.
"""
# Code that pattern matches states a trajectory has been through
# and then cluster them into fundamental few pathways.
#
# This is pulled from 'version 6', which expands on reassigning and
# allows for > 9 states.
#

import pickle
from os.path import exists
from shutil import copyfile

import numpy
import pylcs
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from sklearn.metrics import pairwise_distances
from tqdm.auto import tqdm, trange
from timedinput import timedinput

from lpath.extloader import *
from lpath.io import load_file, default_dendrogram_colors

from ._logger import Logger

log = Logger().get_logger(__name__)


def tostr(b):
    """
    Convert a nonstandard string object ``b`` to str with the handling of the
    case where ``b`` is bytes.

    """
    if b is None:
        return None
    elif isinstance(b, bytes):
        return b.decode('utf-8')
    else:
        return str(b)


def condense_string(string: str, n: int) -> str:
    """
    Function that takes in a string and remove any consecutive duplicates.
    Starts from 1, slowly works up to n.

    Parameters
    ----------
    string : str
        Input string.

    n : int
        How many consecutive duplicates to remove.
        0 for none. 1 for just consecutive characters,
        2 for consecutive pairs (e.g., ABABABA --> AB)

    """
    for n_idx in range(1, n+1):
        result = ""
        count = 0
        while count < len(string):
            if (count + n_idx <= len(string) and
                    string[count:count+n_idx] == string[count+n_idx:count+2*n_idx]):
                count += n_idx
            else:
                result += string[count]
                count += 1
        string = result

    return string


def calc_dist(seq1, seq2, dictionary, pbar, condense=0):
    """
    Pattern match and calculate the similarity between two ``state string`` sequences.

    Parameters
    ----------
    seq1 : numpy.ndarray
        First string to be compared.

    seq2 : numpy.ndarray
        Second string to be compared.

    dictionary : dict
        Dictionary mapping ``state_id`` (float/int) to ``state string`` (characters).

    pbar : tqdm.tqdm
        A tqdm.tqdm object for the progress bar.

    condense : int, default: 0
        Set to a positive int to shorten consecutive characters in state strings.

    Returns
    -------
    1 - similarity : float
        Similarity score.

    """
    # Remove all instances of "unknown" state, which is always the last entry in the dictionary.
    seq1 = seq1[seq1 < len(dictionary) - 1]
    seq1_str = "".join(dictionary[x] for x in seq1)
    seq2 = seq2[seq2 < len(dictionary) - 1]
    seq2_str = "".join(dictionary[x] for x in seq2)

    seq1_str = condense_string(seq1_str, condense)
    seq2_str = condense_string(seq2_str, condense)

    km = pylcs.lcs_sequence_length(seq1_str, seq2_str)
    similarity = (2 * km) / (len(seq1_str) + len(seq2_str) - (abs(len(seq1_str) - len(seq2_str)) / 2))

    pbar.update(1)

    return 1 - similarity


def calc_dist_substr(seq1, seq2, dictionary, pbar, condense=0):
    """
    Pattern match and calculate the similarity between two ``state string`` substrings.
    Used when you're comparing segment ids.

    Parameters
    ----------
    seq1 : numpy.ndarray
        First string to be compared.

    seq2 : numpy.ndarray
        Second string to be compared.

    dictionary : dict
        Dictionary mapping ``state_id`` (float/int) to ``state string`` (characters).

    pbar : tqdm.tqdm
        A tqdm.tqdm object for the progress bar.

    condense : int, default: 0
        Set to a positive int to shorten consecutive characters in state strings.

    Returns
    -------
    1 - similarity : float
        Similarity score.

    """
    # Remove all instances of initial/basis states.
    # seq1 = seq1[seq1 > 0]
    seq1_str = "".join(dictionary[x] for x in seq1)
    # seq2 = seq2[seq2 > 0]
    seq2_str = "".join(dictionary[x] for x in seq2)

    seq1_str = condense_string(seq1_str, condense)
    seq2_str = condense_string(seq2_str, condense)

    km = pylcs.lcs_string_length(seq1_str, seq2_str)
    similarity = (2 * km) / (len(seq1_str) + len(seq2_str) - (abs(len(seq1_str) - len(seq2_str)) / 2))

    pbar.update(1)

    return 1 - similarity


def calc_dist_vanilla(seq1, seq2, dictionary, pbar, condense=0):
    """
    Pattern match and calculate the similarity between two ``state string`` sequences.
    This version does not include the reward term for segments of different length.

    Parameters
    ----------
    seq1 : numpy.ndarray
        First string to be compared.

    seq2 : numpy.ndarray
        Second string to be compared.

    dictionary : dict
        Dictionary mapping ``state_id`` (float/int) to ``state string`` (characters).

    pbar : tqdm.tqdm
        A tqdm.tqdm object for the progress bar.

    condense : int, default: 0
        Set to a positive int to shorten consecutive characters in state strings.

    Returns
    -------
    1 - similarity : float
        Similarity score.

    """
    # Remove all instances of "unknown" state, which is always the last entry in the dictionary.
    seq1 = seq1[seq1 < len(dictionary) - 1]
    seq1_str = "".join(dictionary[x] for x in seq1)
    seq2 = seq2[seq2 < len(dictionary) - 1]
    seq2_str = "".join(dictionary[x] for x in seq2)

    seq1_str = condense_string(seq1_str, condense)
    seq2_str = condense_string(seq2_str, condense)

    km = pylcs.lcs_sequence_length(seq1_str, seq2_str)
    similarity = (2 * km) / (len(seq1_str) + len(seq2_str))

    pbar.update(1)

    return 1 - similarity


def calc_dist_substr_vanilla(seq1, seq2, dictionary, pbar, condense=0):
    """
    Pattern match and calculate the similarity between two ``state string`` substrings.
    Used when you're comparing segment ids.
    This version does not include the reward term for segments of different length.

    Parameters
    ----------
    seq1 : numpy.ndarray
        First string to be compared.

    seq2 : numpy.ndarray
        Second string to be compared.

    dictionary : dict
        Dictionary mapping ``state_id`` (float/int) to ``state string`` (characters).

    pbar : tqdm.tqdm
        A tqdm.tqdm object for the progress bar.

    condense : int, default: 0
        Set to a positive int to shorten consecutive characters in state strings.

    Returns
    -------
    1 - similarity : float
        Similarity score.

    """
    # Remove all instances of initial/basis states.
    # seq1 = seq1[seq1 > 0]
    seq1_str = "".join(dictionary[x] for x in seq1)
    # seq2 = seq2[seq2 > 0]
    seq2_str = "".join(dictionary[x] for x in seq2)

    seq1_str = condense_string(seq1_str, condense)
    seq2_str = condense_string(seq2_str, condense)

    km = pylcs.lcs_string_length(seq1_str, seq2_str)
    similarity = (2 * km) / (len(seq1_str) + len(seq2_str))

    pbar.update(1)

    return 1 - similarity


def determine_reassign(reassign_method):
    """
    Argument processing to determine function to reassign trajectories.

    Parameters
    ----------
    reassign_method : str , default: 'reassign_identity'
        String from argument.reassign_identity, straight from argparser.

    Returns
    -------
    reassign : function
        The reassignment function.

    """
    # Dealing with the preset assign_method
    preset_reassign = {
        'reassign_identity': reassign_identity,
        'reassign_statelabel': reassign_statelabel,
        'reassign_custom': reassign_custom,
        'reassign_segid': reassign_segid,
    }

    if reassign_method in preset_reassign.keys():
        reassign = preset_reassign[reassign_method]
    else:
        import sys
        import os
        sys.path.append(os.getcwd())

        reassign = get_object(reassign_method)
        log.info(f'INFO: Replaced reassign() with {reassign_method}')

    return reassign


def determine_metric(match_metric, match_vanilla):
    """
    Argument processing to determine function to reassign trajectories.

    Parameters
    ----------
    match_metric : str , default: 'longest_common_subsequence'
        String from argument.match_metric, straight from argparser.

    match_vanilla : bool, default: False
        Which similarity metric to use. False to use similarity metric
        with reward term.

    Returns
    -------
    metric : function
        The matching  function.

    """
    # Dealing with the preset assign_method
    subsequence_metric = {
        'longest_common_subsequence': calc_dist,
        'longest_common_subsequence_vanilla': calc_dist_vanilla,
    }

    substring_metric = {
        'longest_common_substring': calc_dist_substr,
        'longest_common_substring_vanilla': calc_dist_substr_vanilla,
    }

    preset_metric = {**subsequence_metric, **substring_metric}

    if match_metric in preset_metric.keys():
        metric = preset_metric[match_metric]
    else:
        import sys
        import os
        sys.path.append(os.getcwd())

        metric = get_object(match_metric)
        log.info(f'INFO: Replaced match_metric with {match_metric}')

    # Dealing with cases where you called the non-vanilla versions (`--substring` or `--subsequence`),
    # but also called `--match-reward-off`. The `--match-reward-off` will take priority.
    if match_vanilla is True:
        if metric in subsequence_metric.values():
            metric = calc_dist_vanilla
        elif metric in substring_metric.values():
            metric = calc_dist_substr_vanilla

    return metric


def load_data(file_name):
    """
    Load in the pickle data from ``extract``.

    Parameters
    ----------
    file_name: str
        File name of the pickle object from ``extract``

    Returns
    -------
    data : list
        A list with the data necessary to reassign, as extracted from ``output.pickle``.

    pathways : numpy.ndarray
        An empty array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    """
    data = load_file(file_name, 1)

    npathways = len(data)
    assert npathways, "Pickle object is empty. Are you sure there are transitions?"
    lpathways = max([len(i) for i in data])
    n = len(data[0][0])

    pathways = numpy.zeros((npathways, lpathways, n), dtype=object)
    # This "Pathways" array should be Iter/Seg/State/auxdata_or_pcoord/frame#/weight

    log.debug(f'Loaded pickle object.')

    return data, pathways


def reassign_custom(data, pathways, dictionary, assign_file=None):
    """
    Reclassify/assign frames into different states. This is highly
    specific to the system. If w_assign's definition is sufficient,
    you can proceed with what's made in the previous step
    using ``reassign_identity``.

    In this example, the dictionary maps state idx to its corresponding ``state_string``.
    We suggest using alphabets as states.

    Parameters
    ----------
    data : list
        An array with the data necessary to reassign, as extracted from ``output.pickle``.

    pathways : numpy.ndarray
        An empty array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    dictionary : dict
        An empty dictionary obj for mapping ``state_id`` with ``state string``. The last entry in
        the dictionary should be the "unknown" state.

    assign_file : str, default : None
        A string pointing to the ``assign.h5`` file. Needed as a parameter for all functions,
        but is ignored if it's an MD trajectory.

    Returns
    -------
    dictionary : dict
        A dictionary mapping each ``state_id`` (float/int) with a ``state string`` (character).
        The last entry in the dictionary should be the "unknown" state.

    """
    # Other example for grouping multiple states into one.
    for idx, pathway in enumerate(data):
        # The following shows how you can "merge" multiple states into
        # a single one.
        pathway = numpy.asarray(pathway)
        # Further downsizing... to if pcoord is less than 5
        first_contact = numpy.where(pathway[:, 3] < 5)[0][0]
        for jdx, frame in enumerate(pathway):
            # First copy all columns over
            pathways[idx, jdx] = frame
            # ortho is assigned to state 0
            if frame[2] in [1, 3, 4, 6, 7, 9]:
                frame[2] = 0
            # para is assigned to state 1
            elif frame[2] in [2, 5, 8]:
                frame[2] = 1
            # Unknown state is assigned 2
            if jdx < first_contact:
                frame[2] = 2
            pathways[idx, jdx] = frame

    # Generating a dictionary mapping each state
    dictionary = {0: 'A', 1: 'B', 2: '!'}

    return dictionary


def reassign_statelabel(data, pathways, dictionary, assign_file):
    """
    Use ``assign.h5`` states as is with ``state_labels``. Does not reclassify/assign frames
    into new states.

    In this example, the dictionary maps state idx to its ``state_labels``,
    as defined in the assign.h5. We suggest using alphabets as ``state_labels``
    to allow for more than 9 states.

    Parameters
    ----------
    data : list
        An list with the data necessary to reassign, as extracted from ``output.pickle``.

    pathways : numpy.ndarray
        An empty array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    dictionary : dict
        An empty dictionary obj for mapping ``state_id`` with "state string".

    assign_file : str
        A string pointing to the ``assign.h5`` file. Needed as a parameter, but ignored if it's an MD trajectory.

    Returns
    -------
    dictionary : dict
        A dictionary mapping each ``state_id`` (float/int) with a state string (character).

    """
    for idx, pathway in enumerate(data):
        pathway = numpy.asarray(pathway)
        for jdx, frame in enumerate(pathway):
            pathways[idx, jdx] = frame

    try:
        import h5py
        with h5py.File(assign_file) as f:
            for idx, state in enumerate(f['state_labels'][:]):
                dictionary[idx] = tostr(state)
        dictionary[len(dictionary)] = '!'  # Unknown state
    except ModuleNotFoundError:
        raise ModuleNotFoundError('Could not import h5py. Exiting out.')

    return dictionary


def reassign_segid(data, pathways, dictionary, assign_file=None):
    """
    Use seg ids as state labels.

    Parameters
    ----------
    data : list
        An list with the data necessary to reassign, as extracted from ``output.pickle``.

    pathways : numpy.ndarray
        An empty array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    dictionary : dict
        An empty dictionary obj for mapping ``state_id`` with ``state string``.

    assign_file : str
        A string pointing to the ``assign.h5`` file. Needed as a parameter, but ignored if it's an MD trajectory.

    Returns
    -------
    dictionary : dict
        A dictionary mapping each ``state_id`` (float/int) with a `state string` (character).

    """
    for idx, pathway in enumerate(data):
        pathway = numpy.asarray(pathway)
        for jdx, frame in enumerate(pathway):
            pathways[idx, jdx] = frame  # Copy everything...
            pathways[idx, jdx, 2] = frame[1]  # Replace states with seg_id

    n_states = int(max([seg[2] for traj in pathways for seg in traj])) + 1
    for idx in range(n_states):
        dictionary[idx] = chr(idx + 65)  # Map seg_id to a unique character

    dictionary[n_states] = '!'  # Unknown state

    return dictionary


def reassign_identity(data, pathways, dictionary, assign_file=None):
    """
    Use assign.h5 states as is. Does not attempt to map assignment
    to ``state_labels`` from assign.h5.

    Parameters
    ----------
    data : list
        An list with the data necessary to reassign, as extracted from ``output.pickle``.

    pathways : numpy.ndarray
        An empty array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    dictionary : dict
        An empty dictionary obj for mapping ``state_id`` with ``state string``.

    assign_file : str
        A string pointing to the ``assign.h5`` file. Needed as a parameter, but ignored if it's an MD trajectory.

    Returns
    -------
    dictionary : dict
        A dictionary mapping each ``state_id`` (float/int) with a `state string` (character).

    """
    for idx, pathway in enumerate(data):
        pathway = numpy.asarray(pathway)
        for jdx, frame in enumerate(pathway):
            pathways[idx, jdx] = frame

    n_states = int(max([seg[2] for traj in pathways for seg in traj])) + 1
    for idx in range(n_states):
        dictionary[idx] = str(idx)

    dictionary[n_states] = '!'  # Unknown state

    return dictionary


def process_shorter_traj(pathways, dictionary, threshold_length, remove_ends):
    """
    Assigns a non-state to pathways which are shorter than
    the max length.

    Parameters
    ----------
    pathways : numpy.ndarray or list
        An array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/frame#/weight.

    dictionary: dict
        Maps each state_id to a corresponding string.

    threshold_length: int or float, default: 0
        A parameter such that trajectories < threshold_length are excluded from pattern matching.

    remove_ends: bool, default: False
        If True, remove the first and last frames (source and target frames).

    """
    del_list = []
    empty_row = numpy.zeros(len(pathways[0, 0]))
    empty_row[2] = len(dictionary) - 1
    for idx, pathway in enumerate(pathways):
        count = 0
        for frame in pathway:
            if frame[0] == 0:  # If no iteration number (i.e., a dummy frame)
                frame[2] = len(dictionary) - 1  # Mark with the last entry
            else:
                count += 1
        if count < threshold_length:
            del_list.append(idx)
        if remove_ends:
            # print(f'{idx} {count} : {len(pathways[idx])}')
            # if count > len(pathways[idx]):
            #     print(pathways[idx])
            pathways[idx, 0] = empty_row
            pathways[idx, count-1] = empty_row

    if remove_ends:
        # Delete in one go.
        pathways = numpy.delete(pathways, [0, -1], axis=1)

    if len(del_list) > 0:
        pathways = numpy.delete(pathways, del_list, axis=0)
        log.debug(f'Indices of pathways removed: {del_list}.')
        log.info(f'Removed {len(del_list)} trajectories of length < {threshold_length} frames.')

    return pathways


def gen_dist_matrix(pathways, dictionary, file_name='succ_traj/distmat.npy', remake=True,
                    metric=calc_dist, condense=0, n_jobs=None):
    """
    Generate the path_string to path_string similarity distance matrix.

    Parameters
    ----------
    pathways : numpy.ndarray
        An array with all the sequences to be compared.

    dictionary : dict
        A dictionary to map pathways states to characters.

    file_name : str, default : 'distmat.npy'
        The file to output the distance matrix.

    remake : bool, default : True
        Indicates whether to remake distance matrix or not.

    metric : bool, default : calc_dist
        Metric function to use.

    condense : int, default: 0
        Set to a positive int to shorten consecutive characters in state strings.

    n_jobs : int, default : None
        Number of jobs to run for the pairwise_distances() calculation. The default issues one job.

    Returns
    -------
    distmat : numpy.ndarray
        A condensed form of the distance matrix (Upper Triangle).

    weights : ndarray
        An array of the weights of each successful pathway (as taken from the last frame).

    """
    weights = []
    path_strings = []
    if metric:
        for pathway in pathways:
            # remove weights for "non-existent" iters (unknown state)
            nonzero = pathway[pathway[:, 2] < len(dictionary) - 1]
            weights.append(nonzero[-1][-1])
            # Create path_strings
            path_strings.append(pathway[:, 2])
    else:
        for pathway in pathways:
            weights.append(pathway[-1][-1])
            # Create path_strings
            path_strings.append(pathway[:, 2])

    weights = numpy.asarray(weights)

    if not exists(file_name) or remake is True:
        log.debug(f'Proceeding to calculate distance matrix.')
        pbar = tqdm(total=int((len(path_strings) * (len(path_strings) - 1))), leave=False)
        # sklearn.metrics.pairwise_distances() brute-forces it by running all n^2 calculations!
        distmat = pairwise_distances(
            X=path_strings, metric=lambda x, y: metric(x, y, dictionary, pbar, condense), n_jobs=n_jobs,
        )
        numpy.save(file_name, distmat)
    else:
        distmat = numpy.load(file_name)
        log.debug(f'Loaded precalculated distance matrix.')

    distmat = squareform(distmat, checks=False)

    return distmat, weights


def calc_linkage(distmat):
    """
    Given distance matrix, calculate and return the linkage.

    """
    return sch.linkage(distmat, method='ward')


def visualize(z, threshold, out_path="plots", show_fig=True, mpl_colors=None, ax=None):
    """
    Visualize the Dendrogram to determine hyper-parameters (n-clusters).
    Theoretically done only once to check.

    Returns
    -------
    plt.gca() : matplotlib.Axes
        A matplotlib.Axes object, which should be the axes which is used to plot the dendrogram.

    """
    # Clean slate.
    if ax is None:
        log.debug('Clearing current plot axis for dendrogram.')
        plt.cla()

    # Plot dendrogram
    try:
        sch.set_link_color_palette(mpl_colors[:-1])
        with plt.rc_context({'lines.linewidth': 2}):
            sch.dendrogram(z, no_labels=True, color_threshold=threshold, above_threshold_color=mpl_colors[-1], ax=ax)
    except RecursionError as e:
        # Catch cases where are too many branches in the dendrogram for default recursion to work.
        import sys

        sys.setrecursionlimit(100000)
        log.warning(e)
        log.warning(f'Dendrogram too complex to plot with default settings. Upping the recursion limit.')
        with plt.rc_context({'lines.linewidth': 2}):
            sch.dendrogram(z, no_labels=True, color_threshold=threshold, above_threshold_color=mpl_colors[-1], ax=ax)

    plt.axhline(y=threshold, c="k")
    plt.ylabel("distance")
    plt.xlabel("pathways")
    plt.savefig(f"{out_path}/dendrogram.pdf")
    if show_fig:
        plt.show()

    return plt.gca()


def hcluster(z, n_clusters):
    """
    Scikit-learn Hierarchical Clustering of the different pathways.

    """
    # (Hyper Parameter t=number of cluster)
    cluster_labels = sch.fcluster(z, t=n_clusters, criterion="maxclust")

    return cluster_labels - 1


def determine_clusters(cluster_labels, clusters=None):
    """
    Determine how many clusters to output.

    Parameters
    ----------
    cluster_labels : numpy.ndarray
        An array with cluster assignments for each pathway.

    clusters : list or None
        Straight from the argparser.

    Returns
    -------
    clusters : list
        A list of clusters to output.

    """
    if clusters is None:
        clusters = list(range(max(cluster_labels) + 1))
    elif not isinstance(clusters, list):
        try:
            list(clusters)
        except TypeError:
            raise TypeError(
                "Provided cluster numbers don't work. Provide a list desired of cluster numbers or 'None' to output "
                "all clusters."
            )

    return clusters


def export_pickle(pathways, output_path):
    """
    Option to output the reassigned pickle object.

    Parameters
    ----------
    pathways : numpy.ndarray
        A reassigned pathway object

    output_path : str
        Path to output pickle object.

    """
    with open(output_path, 'wb') as f:
        pickle.dump(pathways, f)


def select_rep(data_arr, weights, cluster_labels, icluster):
    """
    Small function to determine representative array/weight

    Parameters
    ----------
    data_arr : numpy.ndarray
        The array with all the pathways.

    weights : numpy.ndarray
        Weight information of the pathways.

    cluster_labels : numpy.ndarray
        An array with cluster assignments for each pathway.

    icluster : int
        Index of cluster to look at.

    Returns
    -------
    data_cl : list
        A list of pathways from icluster.

    rep_weight : float
        The weight of the representative structure of icluster.

    """
    selected_cluster = cluster_labels == icluster
    cluster_arr = numpy.array(data_arr, dtype=object)[selected_cluster]
    data_cl_trim = [pathway[pathway[:, 0] != 0] for pathway in cluster_arr]
    weights_cl = weights[selected_cluster]
    rep_weight = data_cl_trim[numpy.argmax(weights_cl)][-1]

    return data_cl_trim, rep_weight


def export_std_files(data_arr, weights, cluster_labels, clusters=None, out_dir="succ_traj"):
    """
    Export data for standard simulations.

    Parameters
    ----------
    data_arr : numpy.ndarray
        The array with all the pathways.

    weights : numpy.ndarray
        Weight information of the pathways.

    cluster_labels : numpy.ndarray
        An array with cluster assignments for each pathway.

    clusters : list or None
        A list of clusters to output, straight from the argparser.

    out_dir : str
        Directory to output files.

    """
    clusters = determine_clusters(cluster_labels, clusters)

    representative_file = f'{out_dir}/representative_segments.txt'
    representative_list = []

    for icluster in clusters:
        trace_out_list = []
        data_cl, rep_weight = select_rep(data_arr, weights, cluster_labels, icluster)
        log.info(f'cluster {icluster} representative weight: {rep_weight}')
        representative_list.append(f'{rep_weight}\n')

        for idx, item in enumerate(data_cl):
            trace_out_list.append(list(numpy.array(item)[:, :2]))

    with open(representative_file, 'w') as f:
        f.writelines(representative_list)


def export_we_files(data_arr, weights, cluster_labels, clusters, file_pattern="west_succ_c{}.h5",
                    out_dir="succ_traj", west_name='west.h5'):
    """
    Export each group of successful trajectories into independent west.h5 file.

    Parameters
    ----------
    data_arr : numpy.ndarray
        The array with all the pathways.

    weights : numpy.ndarray
        Weight information of the pathways.

    cluster_labels : numpy.ndarray
        An array with cluster assignments for each pathway.

    clusters : list or None
        A list of clusters to output.

    file_pattern : str
        String pattern of how files should be outputted.

    out_dir : str
        Directory to output files.

    west_name : str
        Name of west.h5 file to use as base.

    """
    try:
        import h5py
    except ModuleNotFoundError:
        raise ModuleNotFoundError('Could not import h5py. Exiting out.')

    clusters = determine_clusters(cluster_labels, clusters)

    representative_file = f'{out_dir}' + '/representative_segments.txt'
    representative_list = []
    for icluster in clusters:
        new_file = f'{out_dir}/' + file_pattern.format(str(icluster))

        if not exists(new_file):
            copyfile(west_name, new_file)

        first_iter = 1
        with h5py.File(west_name, "r") as h5_file:
            last_iter = len(h5_file['iterations'])

        # Identify constituents of a cluster to output.
        trace_out_list = []
        data_cl, rep_weight = select_rep(data_arr, weights, cluster_labels, icluster)

        log.info(f'cluster {icluster} representative weight: {rep_weight}')
        representative_list.append(f'{rep_weight}\n')

        for idx, item in enumerate(data_cl):
            trace_out_list.append(list(numpy.array(item)[:, :2]))

        # tqdm load bar, working backwards
        tqdm_iter = trange(last_iter, first_iter - 1, -1, desc=f'c{icluster} iterations')

        exclusive_set = {tuple(pair) for ilist in trace_out_list for pair in ilist}
        with h5py.File(new_file, "r+") as h5file:
            for n_iter in tqdm_iter:
                for n_seg in trange(len(h5file[f'iterations/iter_{n_iter:>08}/seg_index']),
                                    desc=f'iter {n_iter} segments', leave=False, delay=1):
                    if (n_iter, n_seg) not in exclusive_set:
                        h5file[f"iterations/iter_{n_iter:>08}/seg_index"]["weight", n_seg] = 0

    with open(representative_file, 'w') as f:
        f.writelines(representative_list)


def determine_rerun(z, out_path='plots', mpl_colors=default_dendrogram_colors, ax=None, timeout=None):
    """
    Asks if you want to regenerate the dendrogram.

    Parameters
    ----------
    z : numpy.ndarray
        A numpy.ndarray from sch.linkage.

    out_path : str, default: 'plots'
        Path to output plots.

    mpl_colors : list or default_dendrogram_colors
        A list of colors for coloring the dendrogram.

    ax : matplotlib.Axes, Default: None
        Matplotlib.Axes object to be inherited.

    timeout : int, default: 30
        Input timeout in seconds.

    """
    if timeout is None:
        timeout = 30

    while True:
        try:
            ans = timedinput('Do you want to regenerate the graph with a new threshold (y/[n])?\n',
                             timeout=timeout, default='N')
            if ans == 'y' or ans == 'Y':
                ans2 = timedinput('What new threshold would you like?\n', timeout=15, default=0.5)
                try:
                    ax = visualize(z, out_path=out_path, threshold=float(ans2),
                                   show_fig=True, mpl_colors=mpl_colors, ax=ax)
                    return ax
                except ValueError:
                    determine_rerun(z, out_path=out_path, mpl_colors=mpl_colors, ax=ax, timeout=timeout)
            elif ans == 'n' or ans == 'N' or ans == '':
                return None
            else:
                log.warning("Invalid input.\n")
        except KeyboardInterrupt:
            sys.exit(0)


def ask_number_clusters(num_clusters=None, timeout=None):
    """
    Asks how many clusters you want to separate the trajectories into.

    """
    if timeout is None:
        timeout = 15

    if not num_clusters:
        while True:
            try:
                ans = timedinput('How many clusters would you like to separate the pathways into?\n',
                                 timeout=timeout, default=2)
                try:
                    ans = int(ans)
                    return ans
                except ValueError:
                    log.warning("Invalid input.\n")
            except KeyboardInterrupt:
                sys.exit(0)
    else:
        return num_clusters


def report_statistics(n_clusters, cluster_labels, weights, segid_status=False):
    """
    Report statistics about the final clusters.

    Parameters
    ----------
    n_clusters : int
        Number of clusters.

    cluster_labels : numpy.ndarray
        An array mapping pathways to cluster

    weights : numpy.ndarray
        Weight information

    segid_status : bool, default: False
        Status of whether we're using seg_ids or not.

    """
    # Initialize the dictionary with 0 weight. 0-based for cl.
    final_dictionary = dict()
    counts = dict()
    uniques = dict()
    for j in range(n_clusters):
        final_dictionary[j] = 0
        counts[j] = 0

    for (cl, weight) in zip(cluster_labels, weights):
        final_dictionary[cl] += weight
        counts[cl] += 1

    if segid_status:
        for cluster, unique_count in zip(*numpy.unique(cluster_labels, return_counts=True)):
            uniques[cluster] = unique_count
    else:
        for cl in range(n_clusters):
            uniques[cl] = 'N/A'

    report = '\n'
    report += f'===LPATH Pattern Matching Statistics===\n'
    report += f'   Total Number of clusters: {n_clusters}\n'
    for (key, val) in final_dictionary.items():
        report += f'   Weight/count/unique count of cluster {key}: {val:.8e} / {counts[key]} / {uniques[key]}\n'
    log.info(report)


def main(arguments):
    """
    Main function that executes the whole `match` step.

    Pathways are processed in the following order:
        1. Assign extra "padding" frames as unknown state for segments too short.
        2. Remove pathways that are too short (if specified).
        3. Remove end frames (source and target states).
        4. During pairwise calculation, condense the frames during comparison and remove frames in "unknown" state.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    # Dealing with the preset assign_method
    reassign = determine_reassign(arguments.reassign_method)
    metric = determine_metric(arguments.match_metric, arguments.match_vanilla)

    # Prepping the data + Calculating the distance matrix
    data, pathways = load_data(arguments.extract_output)

    dictionary = {}
    # Reassignment... (or not)
    dictionary = reassign(data, pathways, dictionary, arguments.assign_name)  # system-specific reassignment of states

    if len(dictionary) < 3:
        log.warning(f'Only {len(dictionary)} states defined, including the "unknown" state. '
                    'This will likely produce bad clustering results and you should considering reassigning to more '
                    'intermediate states using a modified ``--reassign-method``.')

    log.debug(f'Completed reassignment.')

    # Cleanup
    test_obj = process_shorter_traj(pathways, dictionary,
                                    arguments.exclude_short, arguments.remove_ends)
    log.debug(f'Cleaned up trajectories.')

    dist_matrix, weights = gen_dist_matrix(test_obj, dictionary, file_name=arguments.dmatrix_save,
                                           remake=arguments.dmatrix_remake,  # Calculate distance matrix
                                           metric=metric,  # Which metric to use
                                           condense=arguments.condense,  # Whether to condense consecutive state strings
                                           n_jobs=arguments.dmatrix_parallel)  # Number of jobs for pairwise_distance

    log.debug(f'Generated distance matrix.')

    # Visualize the Dendrogram and determine how clusters used to group successful trajectories
    z = calc_linkage(dist_matrix)
    ax = visualize(z, threshold=arguments.dendrogram_threshold, out_path=arguments.out_path,
                   show_fig=arguments.dendrogram_show, mpl_colors=arguments.mpl_colors)
    ax = determine_rerun(z, out_path=arguments.out_path, mpl_colors=arguments.mpl_colors, ax=ax,
                         timeout=arguments.plot_timeout)

    n_clusters = ask_number_clusters(arguments.num_clusters, timeout=arguments.plot_timeout)
    cluster_labels = hcluster(z, n_clusters)

    # Report statistics
    if arguments.stats:
        log.debug('Reporting statistics')
        if arguments.reassign_method == 'reassign_segid':
            segid_status = True
        else:
            segid_status = False
        report_statistics(n_clusters, cluster_labels, weights, segid_status)

    # Output cluster labels and reassigned pickle object
    log.info('Outputting files')
    export_pickle(test_obj, arguments.output_pickle)
    numpy.save(arguments.cl_output, cluster_labels)

    # Following exports each cluster to its own h5 file, all weights of segments not in that group = 0.
    if arguments.export_h5:
        export_we_files(
            test_obj,
            weights,
            cluster_labels,
            clusters=arguments.clusters,
            out_dir=arguments.out_dir,
            file_pattern=arguments.file_pattern,
            west_name=arguments.west_name,
        )
    else:
        export_std_files(
            test_obj,
            weights,
            cluster_labels,
            clusters=arguments.clusters,
            out_dir=arguments.out_dir,
        )
