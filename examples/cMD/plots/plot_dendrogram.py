import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from sys import argv

if len(argv) == 3:
    # Grab number of clusters from argument input.
    n_clusters = argv[1]
    threshold = argv[2]
else:
    # Defaults to 2 clusters and threshold = 2.25.
    n_clusters = 2
    threshold = 2.25

plt.style.use("./default.mplstyle")

plt.figure()

color_list = ["tomato", "dodgerblue", "orchid", "orange", "cyan"]

distmat = np.load("../succ_traj/distmat.npy")
distmat_condensed = squareform(distmat, checks=False)
z = sch.linkage(distmat_condensed, method="ward")
labels         = sch.fcluster(z, t=n_clusters, criterion="maxclust") - 1
labels_str     = [f"cluster #{l}: n={c}\n" for (l,c) in zip(*np.unique(labels, return_counts=True))]
n_clusters     = len(labels_str)
cluster_colors_array = [color_list[l] for l in labels]
link_cols = {}
for i, i12 in enumerate(z[:,:2].astype(int)):
    c1, c2 = (link_cols[x] if x > len(z) else cluster_colors_array[x] for x in i12)
    link_cols[i+1+len(z)] = c1 if c1 == c2 else 'grey'
try:
    with plt.rc_context({'lines.linewidth': 3}):
        sch.dendrogram(z, no_labels=True, color_threshold=threshold, link_color_func=lambda x: link_cols[x], above_threshold_color="grey")
except RecursionError as e:
    import sys
    sys.setrecursionlimit(100000)
    log.warning(e)
    log.warning(f'WARNING: Dendrogram too complex to plot with default settings. Upping the recursion limit.')
    with plt.rc_context({'lines.linewidth': 3}):
        sch.dendrogram(z, no_labels=True, color_threshold=threshold, link_color_func=lambda x: link_cols[x], above_threshold_color="grey")

plt.axhline(y=threshold, c="k", linestyle="--", linewidth=2.5)
plt.ylabel("distance")
plt.xlabel("pathways")
plt.tight_layout()
plt.show()
#plt.savefig("dendrogram_edit.pdf")
