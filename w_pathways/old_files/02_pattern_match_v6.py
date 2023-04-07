# 02_pattern_match_v6.py
#
# Code that pattern matches states a trajectory has been through
# and then cluster them into fundamental few pathways.
#
# Version 6 expands on reassigning and allows for > 9 states.
#
# Written by Anthony Bogetti and Jeremy Leung
# Last modified: Jan 30th, 2023

import numpy
import pickle
import pylcs
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
from tqdm.auto import trange
from shutil import copyfile
from os.path import exists
import h5py


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

def calc_dist(seq1, seq2):
    seq1 = seq1[seq1 > -1]
    seq1_str = "".join(dictionary[x] for x in seq1)
    seq2 = seq2[seq2 > -1]
    seq2_str = "".join(dictionary[x] for x in seq2)

    len_seq1 = len(seq1_str)
    len_seq2 = len(seq2_str)

    lcsstr = pylcs.lcs_sequence_length(seq1_str, seq2_str)
    km = int(lcsstr)
    similarity = (2 * km) / (int(len_seq1) + int(len_seq2))

    return 1 - similarity


def load_data(file_name, n=0):
    """
    Load in the pickle data from step 1.

    Parameters
    ==========

    file_name: str
        File name of the pickle object from step 1

    n : int
        The number of extra datasets you would also like to include for 
        later use. The final "pathways" object created has
        iteration/segment/state_label/weights by default.
    """
    with open(file_name, "rb") as f:
        data = pickle.load(f)

    npathways = len(data)
    lpathways = max([len(i) for i in data])

    pathways = numpy.zeros((npathways, lpathways, n + 4), dtype=object)
    # This "Pathways" array should be Iter/Seg/State/auxdata_or_pcoord/weight

    return data, pathways


def reassign(data, pathways, assign_file):
    """
    Reclassify/assign frames into different states. This is highly
    specific to the system. If w_assign's definition is already
    ok, you can proceed with what's made in the previous step.

    In this example, the dictionary maps state idx to its statelabels,
    as defined in the assign.h5. I suggest using alphabets as states.
    """
    for idx, val in enumerate(data):
        flipped_val = numpy.asarray(val)[::-1]
        for idx2, val2 in enumerate(flipped_val):
            pathways[idx, idx2] = val2

    # Other example for grouping multiple states into one.
    # for idx, val in enumerate(data):
    #     flipped_val = numpy.asarray(val, dtype=object)[::-1]
    #     for idx2, val2 in enumerate(flipped_val):
    #         val2[2] = map_dict[val2[2]]
    #         pathways[idx, idx2] = val2
    #         print(val2)
    # 
    #     # The following shows how you can "merge" multiple states into 
    #     # a single one.
    #     flipped_val = numpy.asarray(val)[::-1]
    #     first_contact = numpy.where(flipped_val[:, 3] < 5)[0][0]
    #     for idx2, val2 in enumerate(flipped_val):
    #         # ortho is assigned to state 1
    #         if val2[2] in [1, 3, 4, 6, 7, 9]:
    #             val2[2] = 1
    #         # para is assigned to state 2
    #         elif val2[2] in [2, 5, 8]:
    #             val2[2] = 2
    #         if idx2 < first_contact:
    #             val2[2] = -1
    #         pathways[idx, idx2] = val2

    # Generating a dictionary
    dictionary = {}
    with h5py.File(assign_file) as f:
        for idx, val in enumerate(f['state_labels'][:]):
            dictionary[idx] = tostr(val)
    dictionary[len(dictionary)] = '!' # Unknown state

    return dictionary


def reassign_identity(data, pathways, assign_file):
    """
    Use assign.h5 states as is. Does not reclassify/assign frames
    into new states.
    """
    for idx, val in enumerate(data):
        flipped_val = numpy.asarray(val)[::-1]
        for idx2, val2 in enumerate(flipped_val):
            pathways[idx, idx2] = val2

    dictionary = {}
    with h5py.File(assign_file) as f:
        for idx in range(len(f['state_labels'][:])):
            dictionary[idx] = str(idx)
    dictionary[len(dictionary)] = str(idx+1) # Unknown state

    return dictionary


def expand_shorter_traj(pathways):
    """
    Assigns a non-state to pathways which are shorter than
    the max length.
    """
    for pathway in pathways:
        for step in pathway:
            if step[0] == 0:
                step[2] = -1


def gen_dist_matrix(
    pathways, file_name="distmap.npy", out_dir="succ_traj", remake=False
):
    """
    Generate the path_string to path_string similarity distance matrix.
    """
    out_dir = f'{out_dir.rsplit("/",1)[0]}'
    new_name = f"{out_dir}/{file_name}"

    weights = []
    path_strings = []
    for pathway in pathways:
        # weights for non-existent iters
        nonzero = pathway[pathway[:, 2] > -1]
        weights.append(nonzero[-1][-1])
        # Create path_strings
        path_strings.append(pathway[:, 2])

    weights = numpy.asarray(weights)

    if not exists(new_name) or remake is True:
        distmat = pairwise_distances(
            X=path_strings, metric=lambda X, Y: calc_dist(X, Y)
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
    out_dir = f'{out_dir.rsplit("/",1)[0]}'

    distmat_condensed = squareform(distmat, checks=False)

    Z = sch.linkage(distmat_condensed, method="ward")

    sch.dendrogram(Z, no_labels=True, color_threshold=threshold)

    plt.axhline(y=threshold, c="k")
    plt.ylabel("distance")
    plt.xlabel("pathway")
    plt.savefig(f"{out_dir}/dendrogram10.pdf")
    if show:
        plt.show()


def hcluster(distmat, n_clusters):
    """
    Scikit-learn Hierarchical Clustering of the different pathways.
    """
    distmat_condensed = squareform(distmat, checks=False)

    Z = sch.linkage(distmat_condensed, method="ward")

    # (Hyper Parameter t=number of clustesr)
    cluster_labels = sch.fcluster(Z, t=n_clusters, criterion="maxclust")

    return cluster_labels


def export_files(
    data_arr,
    weights,
    cluster_labels,
    clusters=None,
    file_pattern="west_succ_c{}.h5",
    out_dir="succ_traj",
    main_file="west.h5",
    assign_file="assign.h5",
):
    """
    Export each group of successful trajectories into independent west.h5 file.
    """
    if clusters is None:
        clusters = list(range(1, max(cluster_labels) + 1))
    elif isinstance(clusters) is not list:
        try:
            list(clusters)
        except:
            raise TypeError(
                "Provided cluster numbers don't work. Provde a list desired of cluster numbers or 'None' to output all clusters."
            )

    representative_file = f'{out_dir.rsplit("/",1)[0]}' + '/representative_segments.txt'
    representative_list = []
    for icluster in clusters:

        new_file = f'{out_dir.rsplit("/",1)[0]}/' + file_pattern.format(str(icluster))

        if not exists(new_file):
            copyfile(main_file, new_file)

        first_iter = 1
        with h5py.File(assign_file, "r") as afile:
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
        with h5py.File(new_file, "r+") as h5file, h5py.File(assign_file, "r") as afile:
            for n_iter in tqdm_iter:
                for n_seg in range(afile["nsegs"][n_iter - 1]):
                    if (n_iter, n_seg) not in exclusive_set:
                        h5file[f"iterations/iter_{n_iter:>08}/seg_index"]["weight", n_seg] = 0

    with open(representative_file,'w') as f:
        f.writelines(representative_list)


def determine_rerun():
    """
    Asks if you want to regenerate the dendrogram.
    """
    while True:
        ans = input('Do you want to regenerate the graph with a new threshold (y/[n])?\n')
        if ans == 'y' or ans == 'Y':
            ans2 = input('What new threshold would you like?\n')
            visualize(dist_matrix, threshold=float(ans2),show=True)
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


if __name__ == "__main__":
    # Prepping the data + Calculating the distance matrix
    data, pathways = load_data("./succ_traj/output.pickle", n=2) # Two extra dimensions for Phi/Psi, Iter/seg/state_id/weight are assumed to be present

    # Reassignment... (or not) Make sure `dictionary` is declared globally since calc_distances() requires it.
    dictionary = reassign(data, pathways, assign_file="ANALYSIS/ALL/assign.h5") # system-specific reassignment of states
    # mapping states > 9 states to their statelabels instead. 
    # dictionary = reassign_identity(data, pathways, assign_file="ANALYSIS/ALL/assign.h5") # Use assign.h5 states as is...

    expand_shorter_traj(data) # Necessary if pathways are of variable length
    dist_matrix, weights = gen_dist_matrix(pathways, out_dir="succ_traj", remake=True) # Calculate distance matrix between pathways

    # Visualize the Dendrogram and determine how clusters used to group successful trajectories
    visualize(dist_matrix, threshold=0.5, show=True) # Visualize
    determine_rerun()
    ncluster = ask_number_cluster()
    cluster_labels = hcluster(dist_matrix, ncluster)

    numpy.save('cluster_labels.npy', cluster_labels)

    # Following exports each cluster to its own h5 file, all weights of segments not in that group = 0.
    export_files(
       data,
       weights,
       cluster_labels,
       clusters=None,
       out_dir="succ_traj",
       file_pattern="west_succ_c{}.h5",
       main_file="multi.h5",
       assign_file="ANALYSIS/ALL/assign.h5",
    )
