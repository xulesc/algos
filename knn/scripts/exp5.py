# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 03:57:27 2014

@author: xule
"""

def read_data(fname):
    fvs = []
    for l in open(fname, 'r'):
        d = l.split(' ')
        fvs.append([float(x) for x in d[2:]])
    return fvs

tp_fvs = read_data('data/diabetic_data.csv.clust.tp')
fp_fvs = read_data('data/diabetic_data.csv.clust.fp')
fn_fvs = read_data('data/diabetic_data.csv.clust.fn')

from sklearn.manifold import MDS
import numpy as np
tp_mds = MDS(n_components=2)
tp_fvs_trans = tp_mds.fit_transform(np.array(tp_fvs))

fp_mds = MDS(n_components=2)
fp_fvs_trans = fp_mds.fit_transform(np.array(fp_fvs))

fn_mds = MDS(n_components=2)
fn_fvs_trans = fn_mds.fit_transform(np.array(fn_fvs))

import matplotlib.pyplot as plt
d1 = np.transpose(tp_fvs_trans)
d2 = np.transpose(fp_fvs_trans)
d3 = np.transpose(fn_fvs_trans)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(d2[0], d2[1], s=10, c='r', marker="o", label='false positives')
ax1.scatter(d3[0], d3[1], s=10, c='y', marker="o", label='false negatives')
ax1.scatter(d1[0], d1[1], s=10, c='g', marker="o", label='true positives')
plt.legend(loc='upper left');
plt.show()

import matplotlib.gridspec as gridspec
d1 = np.transpose(np.array(tp_fvs))
d2 = np.transpose(np.array(fp_fvs))
d3 = np.transpose(np.array(fn_fvs))
fig = plt.figure()
sz = 10#len(d1)
grid = gridspec.GridSpec(sz, sz)
for i in xrange(sz):
    x = i + 1
    for j in xrange(sz):
        y = j + 1    
        ax = plt.subplot(grid[i, j])
        ax.scatter(d2[i], d2[j], s=10, c='r', marker="o", label='false positives')
        ax.scatter(d3[i], d3[j], s=10, c='y', marker="o", label='false negatives')
        ax.scatter(d1[i], d1[j], s=10, c='g', marker="o", label='true positives')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.7, 0.3, ('%d-%d' %(i, j)).lstrip('0'), size=8, horizontalalignment='right')
        fig.add_subplot(ax)

all_axes = fig.get_axes()
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(d1[7], d1[6], s=10, c='g', marker="o", label='true positives')
ax1.scatter(d2[7], d2[6], s=10, c='r', marker="o", label='false positives')
ax1.scatter(d3[7], d3[6], s=10, c='y', marker="o", label='false negatives')
plt.legend(loc='upper left');
plt.show()

##
