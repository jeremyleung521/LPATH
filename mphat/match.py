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
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from sklearn.metrics import pairwise_distances
from tqdm.auto import trange

from mphat.extloader import *

log = logging.getLogger(__name__)


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


def calc_dist(seq1, seq2, dictionary):
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

    Returns
    -------
    1 - similarity : float
        Similarity score.

    """
    # Remove all instances of "unknown" state, which is always the last entry in the dictionary.
    seq1 = seq1[seq1 < len(dictionary)]
    seq1_str = "".join(dictionary[x] for x in seq1)
    seq2 = seq2[seq2 < len(dictionary)]
    seq2_str = "".join(dictionary[x] for x in seq2)

    km = int(pylcs.lcs_sequence_length(seq1_str, seq2_str))
    similarity = (2 * km) / (int(len(seq1_str) + len(seq2_str)))

    return 1 - similarity


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
    with open(file_name, "rb") as f:
        data = pickle.load(f)

    npathways = len(data)
    assert npathways, "Pickle object is empty. Are you sure there are transitions?"
    lpathways = max([len(i) for i in data])
    n = len(data[0][0])

    pathways = numpy.zeros((npathways, lpathways, n), dtype=object)
    # This "Pathways" array should be Iter/Seg/State/auxdata_or_pcoord/weight

    return data, pathways


def reassign_custom(data, pathways, dictionary, assign_file=None):
    """
    Reclassify/assign frames into different states. This is highly
    specific to the system. If w_assign's definition is sufficient,
    you can proceed with what's made in the previous step
    using ``reassign_identity``.

    In this example, the dictionary maps state idx to its corresponding ``state_string``.
    I suggest using alphabets as states.

    Parameters
    ----------
    data : list
        An array with the data necessary to reassign, as extracted from ``output.pickle``.

    pathways : numpy.ndarray
        An empty array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/weight.

    dictionary : dict
        An empty dictionary obj for mapping ``state_id`` with ``state string``.

    assign_file : str, default : None
        A string pointing to the ``assign.h5`` file. Needed as a parameter for all functions,
        but is ignored if it's an MD trajectory.

    Returns
    -------
    dictionary : dict
        A dictionary mapping each ``state_id`` (float/int) with a ``state string`` (character).

    """
    # Other example for grouping multiple states into one.
    for idx, val in enumerate(data):
        # The following shows how you can "merge" multiple states into
        # a single one.
        flipped_val = numpy.asarray(val)[::-1]
        # Further downsizing... to if pcoord is less than 5
        first_contact = numpy.where(flipped_val[:, 3] < 5)[0][0]
        for idx2, val2 in enumerate(flipped_val):
            # ortho is assigned to state 0
            if val2[2] in [1, 3, 4, 6, 7, 9]:
                val2[2] = 0
            # para is assigned to state 1
            elif val2[2] in [2, 5, 8]:
                val2[2] = 1
            # Unknown state is assigned 2
            if idx2 < first_contact:
                val2[2] = 2
            pathways[idx, idx2] = val2

    # Generating a dictionary mapping each state
    dictionary = {0: 'A', 1: 'B', 2: '!'}

    return dictionary


def reassign_statelabel(data, pathways, dictionary, assign_file):
    """
    Use ``assign.h5`` states as is with ``statelabels``. Does not reclassify/assign frames
    into new states.

    In this example, the dictionary maps state idx to its ``statelabels``,
    as defined in the assign.h5. I suggest using alphabets as ``statelabels``
    to allow for more than 9 states.

    Parameters
    ----------
    data : list
        An list with the data necessary to reassign, as extracted from ``output.pickle``.

    pathways : numpy.ndarray
        An empty array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/weight.

    dictionary : dict
        An empty dictionary obj for mapping ``state_id`` with "state string".

    assign_file : str
        A string pointing to the ``assign.h5`` file. Needed as a parameter, but ignored if it's an MD trajectory.

    Returns
    -------
    dictionary : dict
        A dictionary mapping each ``state_id`` (float/int) with a state string (character).

    """
    for idx, val in enumerate(data):
        flipped_val = numpy.asarray(val)[::-1]
        for idx2, val2 in enumerate(flipped_val):
            pathways[idx, idx2] = val2

    try:
        import h5py
        with h5py.File(assign_file) as f:
            for idx, val in enumerate(f['state_labels'][:]):
                dictionary[idx] = tostr(val)
        dictionary[len(dictionary)] = '!'  # Unknown state
    except ModuleNotFoundError:
        raise ModuleNotFoundError('Could not import h5py. Exiting out.')

    return dictionary


def reassign_identity(data, pathways, dictionary, assign_file=None):
    """
    Use assign.h5 states as is. Does not attempt to map assignment
    to ``statelabels`` from assign.h5.

    Parameters
    ----------
    data : list
        An list with the data necessary to reassign, as extracted from ``output.pickle``.

    pathways : numpy.ndarray
        An empty array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/weight.

    dictionary : dict
        An empty dictionary obj for mapping ``state_id`` with ``state string``.

    assign_file : str
        A string pointing to the ``assign.h5`` file. Needed as a parameter, but ignored if it's an MD trajectory.

    Returns
    -------
    dictionary : dict
        A dictionary mapping each ``state_id`` (float/int) with a `state string` (character).

    """
    for idx, val in enumerate(data):
        flipped_val = numpy.asarray(val)[::-1]
        for idx2, val2 in enumerate(flipped_val):
            pathways[idx, idx2] = val2

    n_states = int(max([seg[2] for traj in pathways for seg in traj]))
    for idx in range(n_states):
        dictionary[idx] = str(idx)

    dictionary[n_states] = '!'  # Unknown state

    return dictionary


def expand_shorter_traj(pathways, dictionary):
    """
    Assigns a non-state to pathways which are shorter than
    the max length.

    Parameters
    ----------
    pathways : numpy.ndarray or list
        An array with shapes for iter_id/seg_id/state_id/pcoord_or_auxdata/weight.

    dictionary: dict
        Maps each state_id to a corresponding string.

    """
    for pathway in pathways:
        for step in pathway:
            if step[0] == 0:  # If no iteration number (i.e., a dummy frame)
                step[2] = len(dictionary)  # Mark with the last entry


def gen_dist_matrix(pathways, dictionary, file_name="distmap.npy", out_dir="succ_traj", remake=False):
    """
    Generate the path_string to path_string similarity distance matrix.

    """
    out_dir = f'{out_dir.rsplit("/", 1)[0]}'
    new_name = f"{out_dir}/{file_name}"

    weights = []
    path_strings = []
    for pathway in pathways:
        # weights for non-existent iters
        nonzero = pathway[pathway[:, 2] < len(dictionary)]
        weights.append(nonzero[-1][-1])
        # Create path_strings
        path_strings.append(pathway[:, 2])

    weights = numpy.asarray(weights)

    if not exists(new_name) or remake is True:
        distmat = pairwise_distances(
            X=path_strings, metric=lambda x, y: calc_dist(x, y, dictionary)
        )
        numpy.save(file_name, distmat)

    else:
        distmat = numpy.load(file_name)

    return distmat, weights


def visualize(distmat, threshold, out_dir="succ_traj", show=True):
    """
    Visualize the Dendrogram to determine hyper-parameters (n-clusters).
    Theoretically done only once to check.

    """
    out_dir = f'{out_dir.rsplit("/", 1)[0]}'

    distmat_condensed = squareform(distmat, checks=False)

    z = sch.linkage(distmat_condensed, method="ward")

    try:
        sch.dendrogram(z, no_labels=True, color_threshold=threshold)
    except RecursionError as e:
        # Catch cases where are too many branches in the dendrogram for default recursion to work.
        import sys
        sys.setrecursionlimit(100000)
        log.warning(e)
        log.warning(f'WARNING: Dendrogram too complex to plot with default settings. Upping the recursion limit.')
        sch.dendrogram(z, no_labels=True, color_threshold=threshold)

    try:
        import matplotlib.pyplot as plt
    except (ModuleNotFoundError, ImportError) as e:
        log.debug(e)
        log.debug(f"Can not import matplotlib.")
        return

    plt.axhline(y=threshold, c="k")
    plt.ylabel("distance")
    plt.xlabel("pathway")
    plt.savefig(f"{out_dir}/dendrogram.pdf")
    if show:
        plt.show()


def hcluster(distmat, n_clusters):
    """
    Scikit-learn Hierarchical Clustering of the different pathways.

    """
    distmat_condensed = squareform(distmat, checks=False)

    z = sch.linkage(distmat_condensed, method="ward")

    # (Hyper Parameter t=number of clustesr)
    cluster_labels = sch.fcluster(z, t=n_clusters, criterion="maxclust")

    return cluster_labels


def export_files(
        data_arr,
        weights,
        cluster_labels,
        clusters=None,
        file_pattern="west_succ_c{}.h5",
        out_dir="succ_traj",
        west_name="west.h5",
        assign_name="assign.h5",
):
    """
    Export each group of successful trajectories into independent west.h5 file.

    """
    # TODO: Add option to output pathways object, with might be different from output.pickle after reassignment.
    try:
        import h5py
    except ModuleNotFoundError:
        raise ModuleNotFoundError('Could not import h5py. Exiting out.')

    if clusters is None:
        clusters = list(range(1, max(cluster_labels) + 1))
    elif not isinstance(clusters, list):
        try:
            list(clusters)
        except TypeError:
            raise TypeError(
                "Provided cluster numbers don't work. Provide a list desired of cluster numbers or 'None' to output \
                all clusters."
            )

    representative_file = f'{out_dir.rsplit("/", 1)[0]}' + '/representative_segments.txt'
    representative_list = []
    for icluster in clusters:

        new_file = f'{out_dir.rsplit("/", 1)[0]}/' + file_pattern.format(str(icluster))

        if not exists(new_file):
            copyfile(west_name, new_file)

        first_iter = 1
        with h5py.File(assign_name, "r") as afile:
            last_iter = len(afile["nsegs"])

        tqdm_iter = trange(last_iter, first_iter - 1, -1, desc="iter")

        trace_out_list = []
        selected_cluster = cluster_labels == icluster
        cluster_arr = numpy.array(data_arr, dtype=object)[selected_cluster]
        data_cl = list(cluster_arr)
        weights_cl = weights[selected_cluster]
        print(data_cl[numpy.argmax(weights_cl)][0])
        representative_list.append(str(data_cl[numpy.argmax(weights_cl)][0]) + '\n')

        for idx, item in enumerate(data_cl):
            trace_out_list.append(list(numpy.array(item)[:, :2]))

        exclusive_set = {tuple(pair) for ilist in trace_out_list for pair in ilist}
        with h5py.File(new_file, "r+") as h5file, h5py.File(assign_name, "r") as afile:
            for n_iter in tqdm_iter:
                for n_seg in range(afile["nsegs"][n_iter - 1]):
                    if (n_iter, n_seg) not in exclusive_set:
                        h5file[f"iterations/iter_{n_iter:>08}/seg_index"]["weight", n_seg] = 0

    with open(representative_file, 'w') as f:
        f.writelines(representative_list)


def determine_rerun(dist_matrix):
    """
    Asks if you want to regenerate the dendrogram.

    Parameters
    ----------
    dist_matrix : numpy.ndarray
        A Numpy array of the distance matrix.
    """
    while True:
        ans = input('Do you want to regenerate the graph with a new threshold (y/[n])?\n')
        if ans == 'y' or ans == 'Y':
            ans2 = input('What new threshold would you like?\n')
            visualize(dist_matrix, threshold=float(ans2), show=True)
        elif ans == 'n' or ans == 'N' or ans == '':
            break
        else:
            input("Invalid input.\n")


def ask_number_cluster():
    """
    Asks how many clusters you want to separate the trajectories into.

    """
    while True:
        ans = input('How many clusters would you like to separate the pathways into?\n')
        try:
            ans = int(ans)
            return ans
        except ValueError:
            print("Invalid input.\n")


def main(arguments):
    """
    Main function that executes the whole `match` step.

    Parameters
    ----------
    arguments : argparse.Namespace
        A Namespace object will all the necessary parameters.

    """
    # Dealing with the preset assign_method
    preset_reassign = {
        'reassign_identity': reassign_identity,
        'reassign_statelabel': reassign_statelabel,
        'reassign_custom': reassign_custom,
    }

    if arguments.reassign_method in preset_reassign.keys():
        reassign = preset_reassign[arguments.reassign_method]
    else:
        import sys
        import os
        sys.path.append(os.getcwd())

        reassign = get_object(arguments.reassign_method)
        log.info(f'INFO: Replaced reassign() with {arguments.reassign_method}')

    # Prepping the data + Calculating the distance matrix
    data, pathways = load_data(arguments.input_pickle)

    dictionary = {}
    # Reassignment... (or not) Make sure `dictionary` is declared globally since calc_distances() requires it.
    dictionary = reassign(data, pathways, dictionary, arguments.assign_name)  # system-specific reassignment of states

    # Cleanup
    expand_shorter_traj(pathways, dictionary)  # Necessary if pathways are of variable length
    dist_matrix, weights = gen_dist_matrix(pathways, dictionary, file_name=arguments.dmatrix_save,
                                           out_dir=arguments.out_dir,
                                           remake=arguments.dmatrix_remake)  # Calculate distance matrix

    # Visualize the Dendrogram and determine how clusters used to group successful trajectories
    visualize(dist_matrix, threshold=arguments.dendrogram_threshold, out_dir=arguments.out_dir,
              show=arguments.dendrogram_show)  # Visualize
    determine_rerun(dist_matrix)
    ncluster = ask_number_cluster()
    cluster_labels = hcluster(dist_matrix, ncluster)

    # Output cluster labels
    numpy.save(arguments.cl_output, cluster_labels)

    # Following exports each cluster to its own h5 file, all weights of segments not in that group = 0.
    if arguments.we and arguments.export_h5:
        export_files(
            data,
            weights,
            cluster_labels,
            clusters=arguments.clusters,
            out_dir=arguments.out_dir,
            file_pattern=arguments.file_pattern,
            west_name=arguments.west_name,
            assign_name=arguments.assign_name,
        )


if __name__ == "__main__":
    """
    For calling `match.py` directly. Note all of the parameters are specified manually here.
    
    """
    import argparse
    args = argparse.Namespace(
        input_pickle='succ_traj/output.pickle',  # Input file name of the pickle from `extract.py`
        west_name='multi.h5',  # Name of input HDF5 file (e.g., west.h5)
        assign_name='ANALYSIS/ALL/assign.h5',  # Name of input assign.h5 file
        dmatrix_remake=True,  # Enable to remake the distance Matrix
        dmatrix_save='distmap.npy',  # If dmatrix_remake is False, load this file instead. Assumed located in {out_dir}.
        dendrogram_threshold=0.5,  # Threshold for the Dendrogram
        dendrogram_show=True,  # Show the Dendrogram using plt.show()
        out_dir='succ_traj',  # Output for the distance Matrix
        cl_output='succ_traj/cluster_labels.npy',  # Output path for cluster labels
        file_pattern='west_succ_c{}.h5',  # Pattern to name cluster files
        clusters=None,  # Cluster index to output... otherwise None --> All
        reassign_method='reassign_identity',  # Reassign method. Could be a module to be loaded.
    )

    main(args)
