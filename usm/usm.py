#!/usr/bin/python

import numpy as np
from scipy.spatial import distance
import zlib

class USM:
  def __init__(self, contact_dist = 6.5, max_coords = 1000, metric = 'euclidean'):
    self.contact_dist = contact_dist
    self.max_coords = max_coords
    self.metric = metric
    
  def _get_cm_file_string(self, X, Y, n_contacts, n_atoms):
    st = '\n'.join(map(lambda x : '%s %s' %(x[0], x[1]), zip(X, Y)))
    return '%d\t# Number of Residues\n%d\t# Number of Contacts\n%s' \
      %(n_atoms, n_contacts, st)
    
  def _compute_distance(self, x, y, xy, yx):
    return max(yx - y, xy - x) / max(x, y)
    
  def _in_memory_compress(self, to_zip):
    return float(len(zlib.compress(to_zip)))
    
  def _make_contact_map(self, coords):
    dist = distance.cdist(coords, coords, self.metric)
    idx1, idx2 = np.where(dist <= self.contact_dist)
    fidx = np.where((idx2-idx1) >= 2)
    return [idx1[fidx], idx2[fidx], len(fidx[0]), len(coords)]

  def get_contact_map(self, coords):
    [X, Y, n_contacts, n_atoms] = self._make_contact_map(
      coords[0:min(len(coords), self.max_coords)])
    return (zip(X, Y), self._get_cm_file_string(X, Y, n_contacts, n_atoms))
  
  def dist(self, cm1, cm2):
    x = self._in_memory_compress(cm1)
    y = self._in_memory_compress(cm2)
    xy = self._in_memory_compress(cm1 + cm2)
    yx = self._in_memory_compress(cm2 + cm1)
    return self._compute_distance(x, y, xy, yx)



