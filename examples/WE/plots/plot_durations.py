import matplotlib.pyplot as plt
import numpy as np
import pickle
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from sys import argv


if len(argv) == 2:
    # Grab number of clusters from argument input.
    n_clusters = int(argv[1])
else:
    # Defaults to 3 clusters.
    n_clusters = 3

colors = ['tomato', 'dodgerblue', 'orchid', 'mediumseagreen', 'darkorange', 'mediumpurple','grey']

cluster_centers = np.load("../centroids.npy")

with open("../succ_traj/reassigned.pickle", "rb") as f:
    data = pickle.load(f)
    print("There are", len(data), "pathways")
    pathways = []
    path_idxs = np.arange(0,len(data))
    for pathway in data:
        pathways.append(pathway)

plt.style.use("./default.mplstyle")

distmat = np.load("../succ_traj/distmat.npy")

distmat_condensed = squareform(distmat, checks=False)

z = sch.linkage(distmat_condensed, method="ward")

labels = sch.fcluster(z, t=n_clusters, criterion="maxclust") - 1

plt.figure()

xs = [0.1 * i for i in range(n_clusters)]

for cidx, cluster in enumerate(range(n_clusters)):

    plt.subplot(-(-n_clusters // 2), 2, cidx+1) # Two a row, doing ceiling division to determine how many rows available.

    path_idxs_c = path_idxs[labels == cluster]

    weights = []
    durations = []

    for idx, pathway in enumerate(pathways):
        if idx in path_idxs_c:

            pathway = np.array(pathway)
            pathway = pathway[pathway[:, -1] > 0]
            duration = int(len(pathway))
            durations.append(duration)
    plt.hist(durations, color=colors[cidx], density=True)
    plt.xlabel("durations")
    plt.ylabel("probability")
    plt.title("class %s"%(cidx+1))

plt.tight_layout()
plt.savefig('duration.pdf', dpi=300)
plt.show()
