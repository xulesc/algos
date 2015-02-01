#!/usr/bin/python

import Bio.PDB
import numpy
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.QCPSuperimposer import QCPSuperimposer
from datetime import datetime
from numpy import dot, sqrt
import sys

RESIDUE_LENGTH = 8
RMSD_THRESH = 0.0000001

def simple_rms(coords1, coords2):
  l = min([coords1.shape[0], coords2.shape[0]])
  diff = coords1[0:l] - coords2[0:l]
  return sqrt(sum(sum(diff * diff)) / l)

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

def super_impose(super_imposer, ref_atoms, sample_atoms):
  len_ref = len(ref_atoms); len_sample = len(sample_atoms)
  ## extract the full structure coordinates
  p1 = numpy.zeros((len_ref, 3))
  p2 = numpy.zeros((len_sample, 3))
  for i in range(0, len_ref):
    p1[i] = ref_atoms[i].get_coord()
  for i in range(0, len_sample):
    p2[i] = sample_atoms[i].get_coord()
  min_rmsd = sys.float_info.max; best_rot = []; best_tran = []
  for idx1 in range(0, len_ref - RESIDUE_LENGTH):
    r1 = ref_atoms[idx1 : idx1 + RESIDUE_LENGTH]
    for idx2 in range(0, len_sample - RESIDUE_LENGTH):
      r2 = sample_atoms[idx2 : idx2 + RESIDUE_LENGTH]
      ## align residues
      [rms, rot, tran] = super_impose_residues(super_imposer, r1, r2)
      # rotate the entire structure 
      p2_on_p1 = dot(p2, rot) + tran
      # calculate goodness of alignment over the entire structure
      full_rms = simple_rms(p1, p2_on_p1)      
      if full_rms < min_rmsd:
        min_rmsd = full_rms
        best_rot = rot
        best_tran = tran
        if min_rmsd < RMSD_THRESH:
          return [min_rmsd, best_rot, best_tran]
  return [min_rmsd, best_rot, best_tran]
  
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
  # Calculate RMSDs
  t0 = datetime.now()
  rmsd1, rot1, tran1 = super_impose(SVDSuperimposer(), ref_atoms, sample_atoms)
  t1 = datetime.now()
  rmsd2, rot2, tran2 = super_impose(QCPSuperimposer(), ref_atoms, sample_atoms)
  t2 = datetime.now()
  (d1, d2) = (t1-t0, t2-t1)
  return [d1.total_seconds()*1000, d2.total_seconds()*1000, rmsd1, rmsd2]
  
def all_to_all_psc(data_dir = "pdb_data"):
  from os import listdir
  from os.path import isfile, join
  filenames = [ f for f in listdir(data_dir) if isfile(join(data_dir,f)) ]
  svd_times = []; qcp_times = []
  for f1 in filenames:
    for f2 in filenames:
      print '%s:%s' %(f1, f2)
      [t1, t2, r1, r2] = align_pair(f1, f2)
      print '%f,%f' %(t1, t2) 
      svd_times.append(t1)
      qcp_times.append(t2)
  svd_t = sum(svd_times); qcp_t = sum(qcp_times)
  print 'Average time using svd %f milliseconds per pair with a total of %f for the experiment' %(svd_t/len(svd_times), svd_t)
  print 'Average time using qcp %f milliseconds per pair with a total of %f for the experiment' %(qcp_t/len(qcp_times), qcp_t)
  
#print align_pair('1aa9.pdb', '1ash.pdb')
all_to_all_psc()
  
  
 
 
 
 
