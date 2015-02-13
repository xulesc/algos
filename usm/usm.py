#!/usr/bin/python

import Bio.PDB.PDBParser
import numpy as np
import os
from scipy.spatial import distance
import zlib

def get_ca_atom_list(model):
  atoms = []
  for chain in model:
    for res in chain:
      try:
        atoms.append(res['CA'])
      except:
        pass
  return atoms
  
def make_contact_map(pdb_parser, pdb_dir, pdomain):
  structure = pdb_parser.get_structure("reference", "%s/%s" %(pdb_dir, pdomain))
  ## get coordinates of CA atoms from first model
  atoms = np.array(map(lambda x : x.coord, get_ca_atom_list(structure[0])))
  ## find inter-residue contacts
  dist = distance.cdist(atoms, atoms, 'euclidean')
  idx1, idx2 = np.where(dist <= 6.5)
  fidx = np.where((idx2-idx1) >= 2)
  return [idx1, idx2, fidx, len(fidx[0]), len(atoms)]

def get_cm_file_string(X, Y, S, n_contacts, n_atoms):
  st = '\n'.join(map(lambda x : '%s %s' %(x[0], x[1]), zip(X[S], Y[S])))
  return '%d\t# Number of Residues\n%d\t# Numer of Contacts at 6.5 Angstroms\n%s\n' %(n_atoms, n_contacts, st)
  
def write_contact_maps(in_dir, out_dir):
  pdb_parser = Bio.PDB.PDBParser(QUIET = True)
  for file in os.listdir(in_dir):
    [X, Y, S, n_contacts, n_atoms] = make_contact_map(pdb_parser, in_dir, file)  
    f = open('contact_maps/%s.cm' %file, 'w')
    f.write(get_cm_file_string(X, Y, S, n_contacts, n_atoms))
    f.close()
    
def computeDistance(k_xy, k_yx, k_x, k_y):
  return max(k_yx-k_y, k_xy-k_x) / max(k_x, k_y)
  
def in_memory_compress(STUFF_TO_GZIP):
  return float(len(zlib.compress(STUFF_TO_GZIP)))
    
def get_contact_map_complexities(in_dir):
  pdb_parser = Bio.PDB.PDBParser(QUIET = True)
  structure_cm_string = {}; structure_self_complexity = {}

  for file in os.listdir(in_dir):
    [X, Y, S, n_contacts, n_atoms] = make_contact_map(pdb_parser, in_dir, file)  
    structure_cm_string[file] = get_cm_file_string(X, Y, S, n_contacts, n_atoms)
    structure_self_complexity[file] = in_memory_compress(structure_cm_string[file]) 
    
  for k1, v1 in structure_cm_string.iteritems():
    x = structure_self_complexity[k1]
    for k2, v2 in structure_cm_string.iteritems():
      y = structure_self_complexity[k2]
      xy = in_memory_compress(structure_cm_string[k1]+structure_cm_string[k2]) 
      yx = in_memory_compress(structure_cm_string[k2]+structure_cm_string[k1]) 
      print '%s - %s : %f' %(k1, k2, computeDistance(xy, yx, x, y))  

if __name__ == '__main__':  
  np.set_printoptions(threshold='nan')
  ## make contact maps
  #write_contact_maps('pdb_data', 'contact_maps')
  ## make complexities
  get_contact_map_complexities('pdb_data')



