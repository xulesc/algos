# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 13:24:05 2014

@author: anuj
"""

####################
import csv
import numpy as np
from sklearn import preprocessing
##
csv_reader = lambda x: csv.reader(open(x, 'rb'), delimiter=',')
base_stats = lambda x: [min(x), max(x), np.mean(x), np.std(x)]
min_max_scaler = preprocessing.MinMaxScaler()
##

data_file='/home/xule/workspace/assignment/data/abalone.data'
header=['sex','length','diameter','height','whole','shucked','viscera','shell','rings']
data=zip(*map(lambda x: map(lambda y: float(y), x), csv_reader(data_file)))
#####################
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
#fig = plt.figure()
#grid = gridspec.GridSpec(7, 7)
#for i in xrange(7):
#    x = i + 1
#    for j in xrange(7):
#        y = j + 1    
#        ax = plt.subplot(grid[i, j])
#        ax.scatter(data[x], data[y], c=data[8])
#        ax.set_xticks([])
#        ax.set_yticks([])
#        ax.text(0.7, 0.3, ('%s-%s' %(header[x], header[y])).lstrip('0'), size=8, horizontalalignment='right')
#        fig.add_subplot(ax)
#
#all_axes = fig.get_axes()
##
#from sklearn.decomposition import TruncatedSVD
#red = TruncatedSVD(2)
from sklearn.manifold import MDS
fv = zip(*data[1:len(data) - 1])
t = data[len(data) - 1]
#fv1 = fv.astype(np.float64)
#similarities = euclidean_distances()
#print np.abs(similarities - similarities.T).max()
# Prints 1.7763568394e-15
#mds.fit(data.astype(np.float64))
#similarities = euclidean_distances(fv, fv)
#print np.abs(similarities - similarities.T).max() < 1e-13
red = MDS(2, max_iter=20, dissimilarity="euclidean")
fv_np = np.array(fv)
fv2 = red.fit_transform(fv_np)
plt.scatter(map(lambda x: x[0], fv2), map(lambda x: x[1], fv2), c=t)
##
#np.savetxt("dbscan_mds_reduce.csv", fv2, delimiter=",")
##
#from sklearn.naive_bayes import GaussianNB
#from sklearn.neighbors import KNeighborsClassifier
#names = ["Nearest Neighbors",  "Naive Bayes"]
#classifiers = [ KNeighborsClassifier(29), GaussianNB()]
#TRG_SIZE = 3200
###
#X = zip(*data[:len(data)-1]); Y = data[len(data)-1]
#X_trg = min_max_scaler.fit_transform(fv[:TRG_SIZE]); Y_trg = t[:TRG_SIZE]
#X_tst = min_max_scaler.fit_transform(fv[TRG_SIZE:]); Y_tst = t[TRG_SIZE:]
#for name, classifier in zip(names, classifiers):
#    classifier.fit(X_trg, Y_trg)
#    print '%s: %f' %(name, classifier.score(X_tst, Y_tst))
###############################################################




##
