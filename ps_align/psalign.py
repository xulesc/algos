#!/usr/bin/python

import Bio.PDB
import munkers
import numpy as np
from scipy.spatial import distance
import sys
from datetime import datetime

def get_ca_atom_list(model):
  atoms = []; reses = []
  for chain in model:
    for res in chain:
      try:
        reses.append(res)
        atoms.append(res['CA'])
      except:
        pass
  return (atoms, reses)
  
class DPAlign:
  def __init__(self, w):
    self.w = w
    self.M = np.zeros(w.shape)
    self.M.fill(np.nan)
    self.path = []
    
## M[i][j] = w[i][j] + max(w[i-1][j-1], max(M[i'][j-1] V i' < i), max(M[i-1][j'] V j' < j))
  def align(self, i, j):
    I = i; J = j
    for i in range(I + 1):
      for j in range(J + 1):
        if i == 0 and j == 0:
          self.M[i][j] = self.w[i][j]
          continue
        if i == 0:
          v1 = 0; 
        else:
          v1 = self.M.T[j][0:i].max()
        if j == 0:
          v2 = 0
        else:
          v2 = self.M[i][0:j].max()
        if i == 0 or j == 0:
          v3 = 0
        else:
          v3 = self.w[i-1][j-1]
        self.M[i][j] = self.w[i][j] + max(v1, v2, v3)
    return self.M[I-1][J-1]
    
  def _thread(self, i, j):
    self.path.append((i,j))
    if i == 0 or j == 0:
      return
      
    x = i; y = j
    v = [self.M[i-1][j], self.M[i][j-1], self.M[i-1][j-1]]
    idx = v.index(min(v))
    if idx == 0:
      x = i - 1; 
    elif idx == 1:
      y = j - 1
    elif idx == 2:
      x = i - 1; y = j - 1
    else:
      print 'error'
    self._thread(x, y)    
    
  def thread(self):
    self.path = []
    (i, j) = np.where(self.M==self.M.max())
    self._thread(i[0], j[0])
    self.path.reverse()
    return self.path
    
## PDB domains
pdomain1 = '1aa9.pdb'; pdomain2 = '1ash.pdb'
## parse protein structure files
pdb_parser = Bio.PDB.PDBParser(QUIET = True)
ref_structure = pdb_parser.get_structure("reference", "../pdb_data/%s" %pdomain1)
sample_structure = pdb_parser.get_structure("sample", "../pdb_data/%s" %pdomain2)
## make a list of CA atoms (in the firt model of the structures) to align
(ref_atoms, ref_reses) = get_ca_atom_list(ref_structure[0])
(sample_atoms, sample_reses) = get_ca_atom_list(sample_structure[0])
##
N=100
ref_atoms = ref_atoms[0:N]; ref_reses = ref_reses[0:N]
sample_atoms = sample_atoms[0:N]; sample_reses = sample_reses[0:N]
##
## get coordinates of CA atoms
ref_coords = np.array(map(lambda x: x.get_coord(), ref_atoms))
sample_coords = np.array(map(lambda x: x.get_coord(), sample_atoms))
###################################################################################################################
## Align using bipartite graph techniques based on description given in Wang, Y.; Makedon, F.; Ford, J.; Heng Huang, 
## "A bipartite graph matching framework for finding correspondences between structural elements in two proteins," 
## Engineering in Medicine and Biology Society, 2004. IEMBS '04. 26th Annual International Conference of the IEEE , 
## vol.2, no., pp.2972,2975, 1-5 Sept. 2004
###################################################################################################################
dist_matrix = distance.cdist(ref_coords, sample_coords, 'euclidean').astype(np.int32)
## Non-sequential superposition using the munkres algorithm 
## maximal weight maximal cardinality matching (suitable for global match)
t0 = datetime.now()
cost_matrix = np.array(munkers.run_munkers(dist_matrix, 0)).reshape(dist_matrix.shape)
non_zero = cost_matrix > 0
edges = np.column_stack(np.where(non_zero))
dif = datetime.now() - t0
print "non-sequential threading: %d, %s, %d (msec)" %(len(edges), dist_matrix.shape, dif.total_seconds() * 1000)
## enable if printing large matrices np.set_printoptions(threshold='nan')
print [('%s_%d' %(ref_reses[i].get_resname(), i), '%s-%d' %(sample_reses[j].get_resname(),j)) for i,j in edges]

## Sequential superposition using DP
siml_matrix = dist_matrix.max() - dist_matrix
dpAlign = DPAlign(siml_matrix)
t0 = datetime.now()
dpAlign.align(siml_matrix.shape[0] - 1, siml_matrix.shape[1] - 1)
path = dpAlign.thread()
dif = datetime.now() - t0
print "sequential threading: %d, %s, %d (msec)" %(len(path), dist_matrix.shape, dif.total_seconds() * 1000)
print [('%s_%d' %(ref_reses[i].get_resname(), i), '%s-%d' %(sample_reses[j].get_resname(),j)) for i,j in path]

##
