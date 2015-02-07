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
    for i in range(I):
      for j in range(J):
        if i == 0 and j == 0:
          self.M[i][j] = self.w[i][j]
          continue
        if i == 0 or j == 0:
          v1 = v2 = v3 = 0
        else:
          v1 = self.w[i-1][j-1]
          v2 = self.M.T[j-1][0:i].max()
          v3 = self.M[i-1][0:j].max()
        self.M[i][j] = self.w[i][j] + max(v1, v2, v3)
    return self.M[I-1][J-1]
    
  def thread(self):
    self.path = []
    i, j = np.where(self.M==self.M.max())
    i = i[::-1][0]; j = j[::-1][0]
    self.path.append((i, j))        
    while True:
      if i == -1 or j == -1:
        break
      h = self.M[i][j-1]; v = self.M[i-1][j]; d = self.M[i-1][j-1]
      if d >= v and d >= h:
        self.path.append((i, j))        
        i -= 1; j -= 1
      elif v >= h and v >= d:
        j -= 1
      elif h >= v and h >= d:
        i -=1
      else:
        print 'error'
    self.path = self.path[0:len(self.path)]
    self.path.reverse()
    return np.array(self.path)
    
## PDB domains
pdomain1 = '1aa9.pdb'; pdomain2 = '1ash.pdb'
## parse protein structure files
pdb_parser = Bio.PDB.PDBParser(QUIET = True)
ref_structure = pdb_parser.get_structure("reference", "../pdb_data/%s" %pdomain1)
sample_structure = pdb_parser.get_structure("sample", "../pdb_data/%s" %pdomain2)
## make a list of CA atoms (in the firt model of the structures) to align
(ref_atoms, ref_reses) = get_ca_atom_list(ref_structure[0])
(sample_atoms, sample_reses) = get_ca_atom_list(sample_structure[0])
## get coordinates of CA atoms
ref_coords = np.array(map(lambda x: x.get_coord(), ref_atoms))
sample_coords = np.array(map(lambda x: x.get_coord(), sample_atoms))
###################################################################################################################
## Align using bipartite graph techniques based on description given in Wang, Y.; Makedon, F.; Ford, J.; Heng Huang, 
## "A bipartite graph matching framework for finding correspondences between structural elements in two proteins," 
## Engineering in Medicine and Biology Society, 2004. IEMBS '04. 26th Annual International Conference of the IEEE , 
## vol.2, no., pp.2972,2975, 1-5 Sept. 2004
## Note we are not rotating anything just finding a global alignment
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
dpAlign.align(dist_matrix.shape[0], dist_matrix.shape[1])
path = dpAlign.thread()
dif = datetime.now() - t0
#print np.array2string(dpAlign.w, max_line_width=np.inf)
#print np.array2string(dpAlign.M.astype(np.int32), max_line_width=np.inf)
print "sequential threading: %d, %s, %d (msec)" %(len(path), dist_matrix.shape, dif.total_seconds() * 1000)
print [('%s_%d' %(ref_reses[i].get_resname(), i), '%s-%d' %(sample_reses[j].get_resname(),j)) for i,j in path]

##
