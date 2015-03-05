#!/usr/bin/python

import Bio.PDB
import numpy as np
import sys
from datetime import datetime
from numpy import sum, sqrt, dot
from Bio.PDB.QCPSuperimposer import QCPSuperimposer
from Bio.PDB.PSC.align import Align
from scipy.spatial import distance

RESIDUE_LENGTH = 20
RMS_THRESH = 0.001

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
  
def _superimpose(p1, p2, aligner, type):
    sup = QCPSuperimposer(); align = Align()
    lp2 = p2; previous_rms = sys.maxint; pedges = []
    rms = -1; rot = []; tran = []; edges = []
    while True:
        edges = aligner(p1, lp2)
        if type == 0:
            s = set(edges) ^ set(pedges)
            if len(pedges) != 0 and len(s) < 1.1 * len(pedges):
                break
            
        if len(edges) == 1:
            l = min(p1.shape[0], p2.shape[0])
            edges = [(x, x) for x in range(l)]
        pedges = edges
        a = np.array(map(lambda x : p1[x[0]], edges));  b = np.array(map(lambda x : p2[x[1]], edges));
        sup.set(a, b)
        sup.run()
        rms = sup.get_rms()
        if type != 0:
            if rms > previous_rms:
                rms = previous_rms
                break
            else:
                previous_rms = rms
        rot, tran = sup.get_rotran()
        lp2 = dot(p2, rot) + tran
    return [rms, rot, tran, None, None, edges]
    
def superimpose(p1, p2, type = 0):
    align = Align()
    if type == 0:
        aligner = align.sequential
    else:
        aligner = align.non_sequential
    t0 = datetime.now()
    [rms, rot, tran, p1_r, p2_r, edges] = _superimpose(ref_coords, sample_coords, aligner, type)
    dif = datetime.now()-t0
    print "%s: %d, %s, %d (msec)" %(type, len(edges), [len(ref_atoms), len(sample_atoms)], dif.total_seconds() * 1000)
    print [('%s_%d' %(ref_reses[i].get_resname(), i), '%s-%d' %(sample_reses[j].get_resname(),j)) for i,j in edges]
    print "rmsd: %f" %rms
    print "rot: %s" %rot
    print "tran: %s" %tran
    
def ssuperimpose(p1, p2):
    superimpose(p1, p2)

def nssuperimpose(p1, p2):
    superimpose(p1, p2, 1)

## PDB domains
pdomain1 = '1aa9.pdb'; pdomain2 = '6xia.pdb'
## parse protein structure files
pdb_parser = Bio.PDB.PDBParser(QUIET = True)
ref_structure = pdb_parser.get_structure("reference", "pdb_test_data/%s" %pdomain1)
sample_structure = pdb_parser.get_structure("sample", "pdb_test_data/%s" %pdomain2)
## make a list of CA atoms (in the firt model of the structures) to align
(ref_atoms, ref_reses) = get_ca_atom_list(ref_structure[0])
(sample_atoms, sample_reses) = get_ca_atom_list(sample_structure[0])
## get coordinates of CA atoms
ref_coords = np.array(map(lambda x: x.get_coord(), ref_atoms))
sample_coords = np.array(map(lambda x: x.get_coord(), sample_atoms))
## sequential
ssuperimpose(ref_coords, sample_coords)
## non-sequential
nssuperimpose(ref_coords, sample_coords)                                                                                                    
