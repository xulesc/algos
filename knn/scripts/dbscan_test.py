# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 12:45:22 2014

@author: anuj
"""

import csv
import numpy as np
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def read_data(reader, X_col = 5, Y_col = 0):
    return map(lambda x: map(lambda y: float(y), x), reader)
#################################
csv_reader = lambda x: csv.reader(open(x, 'rb'), delimiter=' ')
data_file  = '/home/xule/workspace/assignment/dbscan_test.dat'
#################################
data = read_data(csv_reader(data_file))
fv = map(lambda x: x[:len(x) - 2], data)
p = map(lambda x: x[len(x) - 2], data)
t = map(lambda x: x[len(x) - 1], data)
#################################
fv2 = MDS(2, max_iter=20, dissimilarity="euclidean").fit_transform(np.array(fv))
#################################
fig = plt.figure()
grid = gridspec.GridSpec(1, 2)
max_xx =  max(map(lambda x: x[0], fv2))
max_yy =  max(map(lambda x: x[1], fv2))
for Y, x, ttl in [(p, 0, 'Predicted'), (t, 1, 'Known')]:
    ax = plt.subplot(grid[0, x])
    ax.scatter(map(lambda x: x[0], fv2), map(lambda x: x[1], fv2), c=Y)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(ttl)
    fig.add_subplot(ax)
all_axes = fig.get_axes()
#################################
#


##
