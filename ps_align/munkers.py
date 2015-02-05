#!/usr/bin/python

import Bio.PDB
import munkers
import numpy as np
from scipy.spatial import distance

## Test
#matrix = np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]]).astype(np.int32)
#print matrix
#print np.array(munkers.run_munkers(matrix, 1)).reshape(matrix.shape)
#print "done"

def get_ca_atom_list(model):
  atoms = []
  for chain in model:
    for res in chain:
      try:
        atoms.append(res['CA'])
      except:
        pass
  return atoms

## PDB domains
pdomain1 = '1aa9.pdb'
pdomain2 = '1ash.pdb'
## parse protein structure files
pdb_parser = Bio.PDB.PDBParser(QUIET = True)
ref_structure = pdb_parser.get_structure("reference", "../pdb_data/%s" %pdomain1)
sample_structure = pdb_parser.get_structure("sample", "../pdb_data/%s" %pdomain2)
## make a list of CA atoms (in the firt model of the structures) to align
ref_atoms = get_ca_atom_list(ref_structure[0])
sample_atoms = get_ca_atom_list(sample_structure[0])
## get coordinates of CA atoms
ref_coords = np.array(map(lambda x: x.get_coord(), ref_atoms))
sample_coords = np.array(map(lambda x: x.get_coord(), sample_atoms))
## get distance matrix
dist_matrix = distance.cdist(ref_coords, sample_coords, 'euclidean').astype(np.int32)
print 'in matrix shape: (%d, %d)' %(dist_matrix.shape[0], dist_matrix.shape[1])
cost_matrix = np.array(munkers.run_munkers(dist_matrix, 1)).reshape(dist_matrix.shape)
#print 'cost matrix:'
#print cost_matrix
