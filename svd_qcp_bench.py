#!/usr/bin/python

import Bio.PDB
import numpy
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.QCPSuperimposer import QCPSuperimposer
from datetime import datetime
from numpy import dot, sqrt
import sys

def super_impose_residues(super_imposer, residue1, residue2):
  ## get atom data in format usable by super imposer
  l = len(residue1)
  fixed_coord = numpy.zeros((l, 3))
  moving_coord = numpy.zeros((l, 3))
  for i in range(0, l):
    fixed_coord[i] = residue1[i].get_coord()
    moving_coord[i] = residue2[i].get_coord()                                                
  # use the superimposer to find best
  super_imposer.set(fixed_coord, moving_coord)
  super_imposer.run()
  #
  rms = super_imposer.get_rms()
  rot, tran = super_imposer.get_rotran()
  return [rms, rot, tran]

def super_impose(super_imposer, ref_atoms, sample_atoms, residue_length):
  len_ref = len(ref_atoms); len_sample = len(sample_atoms)
  t0 = datetime.now(); comparisons = 0
  for idx1 in range(0, len_ref - residue_length):
    r1 = ref_atoms[idx1 : idx1 + residue_length]
    for idx2 in range(0, len_sample - residue_length):
      comparisons += 1
      r2 = sample_atoms[idx2 : idx2 + residue_length]
      [rms, rot, tran] = super_impose_residues(super_imposer, r1, r2)
  t1 = datetime.now()
  d = t1 - t0
  return d.total_seconds()*1000 / comparisons
  
def get_ca_atom_list(model):
  atoms = []
  for chain in model:
    for res in chain:
      try:
        atoms.append(res['CA'])
      except:
        pass
  return atoms

def align_pair(pdomain1, pdomain2):
  # parse the files
  pdb_parser = Bio.PDB.PDBParser(QUIET = True)
  ref_structure = pdb_parser.get_structure("reference", "pdb_data/%s" %pdomain1)
  sample_structure = pdb_parser.get_structure("samle", "pdb_data/%s" %pdomain2)
  # Make a list of CA atoms (in the firt model of the structures) to align
  ref_atoms = get_ca_atom_list(ref_structure[0])
  sample_atoms = get_ca_atom_list(sample_structure[0])
  # 
  upper = min(100, len(ref_atoms), len(sample_atoms))
  for residue_length in range(8, upper):
    svd_time = super_impose(SVDSuperimposer(), ref_atoms, sample_atoms, residue_length)
    qcp_time = super_impose(QCPSuperimposer(), ref_atoms, sample_atoms, residue_length)
    print '%d %f %f' %(residue_length, svd_time, qcp_time)
    sys.stdout.flush()

##
align_pair('1aa9.pdb', '1ash.pdb')
