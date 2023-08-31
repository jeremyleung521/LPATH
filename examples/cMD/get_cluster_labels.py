import numpy
import h5py
import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from sklearn import datasets, neighbors

pcoord1_train = []
pcoord2_train = []

# load in all phi/psi data from cMD simulations
X = numpy.load("phi_psi.npy")[:,0].reshape(-1,1)
Y = numpy.load("phi_psi.npy")[:,1].reshape(-1,1)

# load in training data from WE simulations
with h5py.File("../WE/multi.h5", "r") as f:
    for i in range(250,301):
        istring = "iterations/iter_"+str(i).zfill(8)+"/pcoord"
        ipcoord1 = f[istring][:,-1,0]
        ipcoord2 = f[istring][:,-1,1]
        pcoord1_train.append(ipcoord1)
        pcoord2_train.append(ipcoord2)

X_train = numpy.concatenate(pcoord1_train).reshape(-1,1)
Y_train = numpy.concatenate(pcoord2_train).reshape(-1,1)

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

# agglomerative clustering on the training data
ag = AgglomerativeClustering(n_clusters=None, 
                             linkage="average", 
                             distance_threshold=75).fit(data_train)

labels_train = ag.labels_

uniq_labels = numpy.unique(labels_train)

# plot training data cluster assignments and make sure it looks okay
for idx, i in enumerate(uniq_labels):
    plt.scatter(data_train[:,0][labels_train==i], data_train[:,1][labels_train==i])
    centroid = data_train[labels_train==i].mean(axis=0)
    plt.scatter(centroid[0], centroid[1], color="black", zorder=3)
    plt.annotate(str(i), (centroid[0], centroid[1]))
plt.xlim(-210,150)
plt.ylim(-180,180)
plt.xlabel("Φ")
plt.ylabel("Ψ")

plt.show()

answer = input("continue? (y/n)")
if answer == "n":
    exit()

plt.clf()

# now train a KNN model on the trianing set
knn = neighbors.KNeighborsClassifier(n_neighbors=5)
knn.fit(data_train, labels_train)

# predict cluster labels for the full cMD dataset
labels = knn.predict(data)

# save state labels for LPATH
numpy.save("states.npy", labels)

uniq_labels = numpy.unique(labels)

centroids = []

# plot the full-assigned cMD dataset
for idx, i in enumerate(uniq_labels):
    plt.scatter(data[:,0][labels==i], data[:,1][labels==i])
    centroid = data[labels==i].mean(axis=0)
    centroids.append(centroid)
    plt.scatter(centroid[0], centroid[1], color="black", zorder=3)
    plt.annotate(str(i), (centroid[0], centroid[1])) 
plt.xlim(-210,150)
plt.ylim(-180,180)

centroid_arr = numpy.concatenate(centroids).reshape(-1,2)
numpy.save("centroids.npy", centroid_arr)
print(centroid_arr)

plt.xlabel("$\phi$")
plt.ylabel("$\psi$")

plt.show()
