import matplotlib.pyplot as plt
import numpy as np
import pickle
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from sys import argv

if len(argv) == 2:
    # Grab number of clusters from argument input.
    n_clusters = argv[1]
else:
    # Defaults to 2 clusters.
    n_clusters = 2

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

for cidx, cluster in enumerate(range(0, n_clusters)):

    path_idxs_c = path_idxs[labels==cluster]

    weights = []

    for idx, pathway in enumerate(pathways):
        if idx in path_idxs_c:

            pathway = np.array(pathway)
            pathway = pathway[pathway[:,0]>0]
            weight = pathway[-1,-1]
            weights.append(weight)
    plt.bar(xs[cidx], np.sum(weights), width=0.05, color=colors[cidx])

plt.xlim(xs[0]-0.1, xs[-1]+0.1)
plt.xticks(ticks=xs, labels=[f"class {i}" for i in range(1, n_clusters+1)], rotation=45)
plt.ylabel("probability")
plt.tight_layout()
plt.savefig('histogram.pdf', dpi=300)
plt.show()
