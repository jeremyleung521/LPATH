import numpy
import h5py
import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from sklearn import datasets, neighbors

plt.style.use("./plots/default.mplstyle")
cmap = get_cmap("Set2")
colors = cmap.colors

pcoord1_train = []
pcoord2_train = []
pcoord1 = []
pcoord2 = []

# load in all data from WESTPA simulation file
with h5py.File("multi.h5", "r") as f:
    for i in range(1,301):
        istring = "iterations/iter_"+str(i).zfill(8)+"/pcoord"
        ipcoord1 = f[istring][:,:,0]
        ipcoord2 = f[istring][:,:,1]
        pcoord1.append(ipcoord1)
        pcoord2.append(ipcoord2)

# load in training data from WESTPA simulation file (final 50 iterations)
with h5py.File("multi.h5", "r") as f:
    for i in range(250,301):
        istring = "iterations/iter_"+str(i).zfill(8)+"/pcoord"
        ipcoord1 = f[istring][:,-1,0]
        ipcoord2 = f[istring][:,-1,1]
        pcoord1_train.append(ipcoord1)
        pcoord2_train.append(ipcoord2)

X_train = numpy.concatenate(pcoord1_train).reshape(-1,1)
Y_train = numpy.concatenate(pcoord2_train).reshape(-1,1)
X = numpy.concatenate(pcoord1).reshape(-1,1)
Y = numpy.concatenate(pcoord2).reshape(-1,1)

# shift periodic boundary from 0-360 to -210-150 
for pidx, val in enumerate(X_train):
    if val > 150:
        X_train[pidx] -= 360

for pidx, val in enumerate(Y_train):
    if val > 180:
        Y_train[pidx] -= 360

for pidx, val in enumerate(X):
    if val > 150:
        X[pidx] -= 360

for pidx, val in enumerate(Y):
    if val > 180:
        Y[pidx] -= 360

data_train = numpy.concatenate((X_train,Y_train), axis=1)
data = numpy.concatenate((X,Y), axis=1)

# perform agglomerative clustering on training data
ag = AgglomerativeClustering(n_clusters=None, 
                             linkage="average", 
                             distance_threshold=75).fit(data_train)

labels_train = ag.labels_

uniq_labels = numpy.unique(labels_train)

# plot training data cluster assignments
for idx, i in enumerate(uniq_labels):
    plt.scatter(data_train[:,0][labels_train==i], data_train[:,1][labels_train==i], c=colors[idx])
    centroid = data_train[labels_train==i].mean(axis=0)
    plt.scatter(centroid[0], centroid[1], color="black", zorder=3)
plt.xlim(-210,150)
plt.ylim(-180,180)
plt.xticks([-100, 0, 100])
plt.yticks([-100, 0, 100])
plt.xlabel("Φ")
plt.ylabel("Ψ")

plt.tight_layout()
#plt.savefig("agg_clusters.png", dpi=300)
plt.show()

answer = input("continue? (y/n)")
if answer == "n":
    exit()

plt.clf()

# train KNN model from training data clusters
knn = neighbors.KNeighborsClassifier(n_neighbors=5)
knn.fit(data_train, labels_train)

# create an auxiliary dataset to save cluster assignments in the multi.h5 file
with h5py.File("multi.h5", "a") as f:
    for i in range(1,301):
        pcoord1 = []
        pcoord2 = []
        istring = "iterations/iter_"+str(i).zfill(8)+"/pcoord"
        ipcoord1 = f[istring][:,:,0]
        n_segs = ipcoord1.shape[0]
        pcoord_len = ipcoord1.shape[1]
        ipcoord2 = f[istring][:,:,1]
        pcoord1.append(ipcoord1)
        pcoord2.append(ipcoord2)
        X = numpy.concatenate(pcoord1).reshape(n_segs*pcoord_len,1)
        Y = numpy.concatenate(pcoord2).reshape(n_segs*pcoord_len,1)
        for pidx, val in enumerate(X):
            if val > 150:
                X[pidx] -= 360
        
        for pidx, val in enumerate(Y):
            if val > 180:
                Y[pidx] -= 360
        idata = numpy.concatenate((X,Y), axis=1)
        dsstring = 'iterations/iter_'+str(i).zfill(8)+'/auxdata/labels'
        dsstring_no_labels = 'iterations/iter_'+str(i).zfill(8)+'/auxdata'
        labels = knn.predict(idata).reshape(n_segs, pcoord_len, 1)
        try:
            f.create_dataset(dsstring, data=labels)
        except Exception:
            f[dsstring][...] = labels

# predict labels for the entire WESTPA dataset
labels = knn.predict(data)

# save the cluster labels
numpy.save("labeled_points.npy", labels)

uniq_labels = numpy.unique(labels)

centroids = []

# plot all cluster assignments
for idx, i in enumerate(uniq_labels):
    plt.scatter(data[:,0][labels==i], data[:,1][labels==i], c=colors[idx])
    centroid = data[labels==i].mean(axis=0)
    centroids.append(centroid)
    plt.scatter(centroid[0], centroid[1], color="black", zorder=3)

centroid_arr = numpy.concatenate(centroids).reshape(-1,2)
numpy.save("centroids.npy", centroid_arr)

plt.xlim(-210,150)
plt.ylim(-180,180)
plt.xticks([-100, 0, 100])
plt.yticks([-100, 0, 100])
plt.xlabel("Φ")
plt.ylabel("Ψ")

plt.tight_layout()

#plt.savefig("knn_clusters.png", dpi=300)
plt.show()
