#!/usr/bin/python

import numpy as np
from scipy.spatial import distance
import zlib

class USM:
  def __init__(self, contact_dist = 6.5):
    self.contact_dist = contact_dist
    
  def _get_cm_file_string(self, X, Y, S, n_contacts, n_atoms):
    st = '\n'.join(map(lambda x : '%s %s' %(x[0], x[1]), zip(X[S], Y[S])))
    return '%d\t# Number of Residues\n%d\t# Number of Contacts\n%s\n' %(n_atoms, n_contacts, st)
    
  def _computeDistance(self, k_xy, k_yx, k_x, k_y):
    return max(k_yx-k_y, k_xy-k_x) / max(k_x, k_y)
    
  def _in_memory_compress(self, to_zip):
    return float(len(zlib.compress(to_zip)))
    
  def _make_contact_map(self, coords):
    dist = distance.cdist(coords, coords, 'euclidean')
    idx1, idx2 = np.where(dist <= self.contact_dist)
    fidx = np.where((idx2-idx1) >= 2)
    return [idx1, idx2, fidx, len(fidx[0]), len(coords)]

  def get_contact_map(self, coords):
    [X, Y, S, n_contacts, n_atoms] = self._make_contact_map(coords)
    return (zip(X[S], Y[S]), self._get_cm_file_string(X, Y, S, n_contacts, n_atoms))
  
  def usm(self, cm1, cm2):
    x = self._in_memory_compress(cm1)
    y = self._in_memory_compress(cm2)
    xy = self._in_memory_compress(cm1 + cm2)
    yx = self._in_memory_compress(cm2 + cm1)
    return self._computeDistance(xy, yx, x, y)



