import matplotlib.pyplot as plt
import numpy as np
import pickle
import scipy.cluster.hierarchy as sch
import networkx as nx
from scipy.spatial.distance import squareform
from itertools import groupby
from sys import argv


cluster_centers = np.load("../centroids.npy")

with open("../succ_traj/reassigned.pickle", "rb") as f:
    data = pickle.load(f)
    print("There are", len(data), "pathways")
    pathways = []
    path_idxs = np.arange(0,len(data))
    for pathway in data:
        pathways.append(pathway)

if len(argv) == 2:
    # Grab number of clusters from argument input.
    n_clusters = int(argv[1])
else:
    # Defaults to 2 clusters.
    n_clusters = 2

plt.style.use("./default.mplstyle")

distmat = np.load("../succ_traj/distmat.npy")

distmat_condensed = squareform(distmat, checks=False)

z = sch.linkage(distmat_condensed, method="ward")

labels = sch.fcluster(z, t=n_clusters, criterion="maxclust") - 1

plt.figure(figsize=(14,8))

for cidx, cluster in enumerate(range(n_clusters)):
    if n_clusters > 4:
        n_col = 3
    else:
        n_col = 2

    plt.subplot(-(-n_clusters // n_col), n_col,
                cluster + 1)  # Doing ceiling division to determine how many rows available.

    plt.title("class "+str(cluster+1))

    path_idxs_c = path_idxs[labels == cluster]

    cluster_counts, state_list, edge_list = [], [], []

    for idx, pathway in enumerate(pathways):
        if idx in path_idxs_c:

            pathway = np.array(pathway)
            pathway = pathway[pathway[:, 0] > 0]
            states = pathway[:, 2]
            weights = pathway[:, -1]
            states_last = states[::25]
            weights_last = weights[::25]

            weighted_counts = np.array([np.sum(weights_last[states_last == i]) for i in range(12)], dtype=float)
            cluster_counts.append(weighted_counts)
            condensed_state = [int(k) for k,g in groupby(states)]
            condensed_state.insert(0, 6)
            condensed_state.insert(len(condensed_state), 7)
            pstates = [2, 3, 1, 5]
            istates = [8, 9, 10, 11]

            if condensed_state not in state_list:
                state_list.append(condensed_state)

                for c in range(len(condensed_state)-1):
                    edge = [condensed_state[c], condensed_state[c+1]]
                    if edge == [pstates[1],pstates[0]]:
                        edge_list.append([pstates[1],istates[0]])
                        edge_list.append([istates[1],pstates[0]])
                    elif edge == [pstates[0],pstates[1]]:
                        edge_list.append([pstates[0],istates[1]])
                        edge_list.append([istates[0],pstates[1]])
                    elif edge == [pstates[2],pstates[3]]:
                        edge_list.append([pstates[2],istates[3]])
                        edge_list.append([istates[2],pstates[3]])
                    elif edge == [pstates[3],pstates[2]]:
                        edge_list.append([pstates[3],istates[2]])
                        edge_list.append([istates[3],pstates[2]])
                    else:
                        edge_list.append(edge)

    sum_counts = np.sum(cluster_counts, axis=0)
    sum_counts = sum_counts*10000
    sum_counts[6] = 500
    sum_counts[7] = 500
    sum_counts[8] = 1
    sum_counts[9] = 1
    sum_counts[10] = 1
    sum_counts[11] = 1

    G = nx.DiGraph()
    G.add_nodes_from([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]) # 8 and 9 are start and end, 10-13 are periodic
    G.add_edges_from(edge_list)
    mapping = {0: "0", 1: "1", 2: "2", 3: "3", 4: "4", 5: "5", 6: "A", 7: "B", 8: "8", 9: "9", 10: "10", 11: "11"}
    if cidx == 0:
        colors = ["#FF6347", "#FF6347", "#FF6347", "#FF6347", "#FF6347", "#FF6347", "#D3D3D3", "#D3D3D3", "#FF6347", "#FF6347", "#FF6347", "#FF6347"]
    else:
        colors = ["#1E90FF", "#1E90FF", "#1E90FF", "#1E90FF", "#1E90FF", "#1E90FF", "#D3D3D3", "#D3D3D3", "#1E90FF", "#1E90FF", "#1E90FF", "#1E90FF"]
    pos = {"0": cluster_centers[0], "1": cluster_centers[1], "2": cluster_centers[2], "3": cluster_centers[3], "4": cluster_centers[4], "5": cluster_centers[5], "A": (-140,80), "B": (80,-50), "8": (-110, 250), "9": (-110, -250), "10": (60, 250), "11": (60, -250)}
    isolated = list(nx.isolates(G))
    new_isolated = []
    if isolated:
        sum_counts = np.delete(sum_counts, isolated)
        colors = np.delete(colors, isolated)
    G.remove_nodes_from(list(nx.isolates(G)))
    G = nx.relabel_nodes(G, mapping)
    nx.draw_networkx_nodes(G, pos=pos, node_color=colors, alpha=0.5, node_size=sum_counts)
    nx.draw_networkx_edges(G, pos=pos, arrows=True, node_size=250, width=1.5)
    nx.draw_networkx_labels(G, pos=pos, font_weight='bold')
    plt.xlim(-200,140)
    plt.ylim(-215,215)
plt.savefig('network.pdf', dpi=300)
plt.show()
