#!/usr/bin/python

import Bio.PDB.PDBParser
from usm import USM
import os
import numpy as np

def get_ca_atom_list(model):
  atoms = []
  for chain in model:
    for res in chain:
      try:
        atoms.append(res['CA'])
      except:
        pass
  return atoms
   
def get_contact_map_complexities(in_dir):
  usm = USM()
  pdb_parser = Bio.PDB.PDBParser(QUIET = True)
  structure_cm_string = {}; 

  for file in os.listdir(in_dir):
    structure = pdb_parser.get_structure("reference", "%s/%s" %(in_dir, file))
    coords = np.array(map(lambda x : x.coord, get_ca_atom_list(structure[0])))  
    structure_cm_string[file] = usm.get_contact_map(coords)[1]
    
  for k1, v1 in structure_cm_string.iteritems():
    for k2, v2 in structure_cm_string.iteritems():
      dist = usm.usm(structure_cm_string[k1], structure_cm_string[k2])
      print '%s - %s : %f' %(k1, k2, dist)  

if __name__ == '__main__':  
  np.set_printoptions(threshold='nan')
  get_contact_map_complexities('pdb_data')
