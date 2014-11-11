# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 12:12:31 2014

@author: anuj
"""

print(__doc__)

from time import time
import numpy as np
import pylab as pl

from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn.preprocessing import normalize
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cross_validation import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn import preprocessing
from sklearn.naive_bayes import GaussianNB
from sklearn.naive_bayes import MultinomialNB
from sklearn import linear_model, datasets
from sklearn.decomposition import PCA

klas_file = '/home/anuj/workspace.custom/assignment/data/dataset_diabetes/diabetic_data.csv.klass'
data_file = '/home/anuj/workspace.custom/assignment/data/dataset_diabetes/diabetic_data.csv.num'

klasses = np.genfromtxt(klas_file, skip_header=1)
n_data = np.genfromtxt(data_file, usecols=range(2, 11), skip_header=1,delimiter=',', missing_values='?')

np.random.seed(42)

data = normalize(n_data[~np.isnan(n_data).any(axis=1)], axis = 0)

n_samples, n_features = [len(data), len(data[0])]
n_digits = len(np.unique(klasses))
labels = klasses[~np.isnan(n_data).any(axis=1)]

sample_size = int(10.0 * n_samples / 100)

print("n_digits: %d, \t n_samples %d, \t n_features %d"
      % (n_digits, n_samples, n_features))


print(79 * '_')
print('% 9s' % 'init'
      '    time  inertia    homo   compl  v-meas     ARI AMI  silhouette')


def bench_k_means(estimator, name, data):
    t0 = time()
    estimator.fit(data)
    print('% 9s   %.2fs    %i   %.3f   %.3f   %.3f   %.3f   %.3f'
          % (name, (time() - t0), estimator.inertia_,
             metrics.homogeneity_score(labels, estimator.labels_),
             metrics.completeness_score(labels, estimator.labels_),
             metrics.v_measure_score(labels, estimator.labels_),
             metrics.adjusted_rand_score(labels, estimator.labels_),
             metrics.adjusted_mutual_info_score(labels,  estimator.labels_)))

bench_k_means(KMeans(init='k-means++', n_clusters=n_digits, n_init=10),
              name="k-means++", data=data)

bench_k_means(KMeans(init='random', n_clusters=n_digits, n_init=10),
              name="random", data=data)

# in this case the seeding of the centers is deterministic, hence we run the
# kmeans algorithm only once with n_init=1
pca = PCA(n_components=n_digits).fit(data)
bench_k_means(KMeans(init=pca.components_, n_clusters=n_digits, n_init=1),
              name="PCA-based",
              data=data)
print(79 * '_')

###############################################################################
# Visualize the results on PCA-reduced data

#reduced_data = PCA(n_components=2).fit_transform(data)
#kmeans = KMeans(init='k-means++', n_clusters=n_digits, n_init=10)
#kmeans.fit(reduced_data)
#
## Step size of the mesh. Decrease to increase the quality of the VQ.
#h = .02     # point in the mesh [x_min, m_max]x[y_min, y_max].
#
## Plot the decision boundary. For that, we will assign a color to each
#x_min, x_max = reduced_data[:, 0].min() + 1, reduced_data[:, 0].max() - 1
#y_min, y_max = reduced_data[:, 1].min() + 1, reduced_data[:, 1].max() - 1
#xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
#
## Obtain labels for each point in mesh. Use last trained model.
#Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])
#
## Put the result into a color plot
#Z = Z.reshape(xx.shape)
#pl.figure(1)
#pl.clf()
#pl.imshow(Z, interpolation='nearest',
#          extent=(xx.min(), xx.max(), yy.min(), yy.max()),
#          cmap=pl.cm.Paired,
#          aspect='auto', origin='lower')
#
#pl.plot(reduced_data[:, 0], reduced_data[:, 1], 'k.', markersize=2)
## Plot the centroids as a white X
#centroids = kmeans.cluster_centers_
#pl.scatter(centroids[:, 0], centroids[:, 1],
#           marker='x', s=169, linewidths=3,
#           color='w', zorder=10)
#pl.title('K-means clustering on the digits dataset (PCA-reduced data)\n'
#         'Centroids are marked with white cross')
#pl.xlim(x_min, x_max)
#pl.ylim(y_min, y_max)
#pl.xticks(())
#pl.yticks(())
#pl.show()

###############################################################################
data_train, data_test, labels_train, labels_test = train_test_split(data, labels, test_size=0.20, random_state=42)
neigh = KNeighborsClassifier(n_neighbors=5)
neigh.fit(data_train, labels_train)
#print neigh.score(data_test, labels_test)
pred = neigh.predict(data_test)
cm = confusion_matrix(labels_test, pred)
print(cm)
pl.matshow(cm)
pl.title('Confusion matrix')
pl.colorbar()
pl.ylabel('True label')
pl.xlabel('Predicted label')
pl.show()
###############################################################################
klas_file = '/home/anuj/workspace.custom/assignment/data/dataset_diabetes/diabetic_data.csv.klass'
data_file = '/home/anuj/workspace.custom/assignment/data/dataset_diabetes/diabetic_data.csv.non_num'
data_file2 = '/home/anuj/workspace.custom/assignment/data/dataset_diabetes/diabetic_data.csv.num'

klasses = np.genfromtxt(klas_file, skip_header=1)
n_data = np.genfromtxt(data_file, skip_header=1,delimiter=',',dtype='|S5')
n_data_num = np.genfromtxt(data_file2, usecols=range(2, 11), skip_header=1,delimiter=',', missing_values='?')
#n_data = n_data[~np.isnan(n_data).any(axis=1)]

exc = np.isnan(n_data_num).any(axis=1)
n_data_num_n = normalize(n_data_num[~exc], axis = 0)
labels = klasses[~exc]
n_data2 = n_data[~exc]
n_data2 = [x[:len(x) - 1] for x in n_data2]

n_data2 = np.transpose(n_data2)
le = preprocessing.LabelEncoder()
n_data3 = [le.fit(d).transform(d) for d in n_data2]

##############
#f_data = np.transpose(n_data3)
f_data = n_data_num_n
#for x in np.transpose(n_data_num_n):
#    f_data.append(x)
#f_data = np.transpose(f_data)
##############

data_train, data_test, labels_train, labels_test = train_test_split(f_data, labels, test_size=0.20, random_state=42)
neigh = KNeighborsClassifier(n_neighbors=1)
#neigh = MultinomialNB()
print('%d:%d\n' %(sum(labels_train),len(labels_train)))
neigh.fit(data_train, labels_train)
#print neigh.score(data_test, labels_test)
pred = neigh.predict(data_test)
print('%d:%d:%d:%d\n' %(sum(labels_test),len(labels_test),sum(pred),len(pred)))
cm = confusion_matrix(labels_test, pred)
print(cm)
pl.matshow(cm)
pl.title('Confusion matrix')
pl.colorbar()
pl.ylabel('True label')
pl.xlabel('Predicted label')
pl.show()
##############
###############################################################################
f = '/home/anuj/workspace.custom/assignment/data/dataset_diabetes/diabetic_data.csv'
d1 = np.genfromtxt(f, delimiter = ',', names=True)
names = d1.dtype.names
d = np.genfromtxt(f, delimiter = ',', dtype='|S5', skip_header=1)
dc = np.transpose(d)
for x, name in zip(dc, names):
    print '%s: %d' %(name, len(np.unique(x)))
##

