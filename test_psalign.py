#!/usr/bin/python

import Bio.PDB
import numpy as np
import sys, os
from datetime import datetime
from numpy import sum, sqrt, dot
from Bio.PDB.QCPSuperimposer import QCPSuperimposer
from Bio.PDB.PSC.align import Align
from Bio.PDB.Polypeptide import PPBuilder
from scipy.spatial import distance
import multiprocessing
import random
from Bio import pairwise2

RESIDUE_LENGTH = 20
RMS_THRESH = 0.01
align = Align()
aligner = [align.sequential, align.non_sequential]
sup = QCPSuperimposer()

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
  
def gentask(p1, p2, type):
    i = 0
    while i < p1.shape[0] - RESIDUE_LENGTH:
        j = 0
        while j < p2.shape[0] - RESIDUE_LENGTH:
            yield [p1,p2,i,j,type]
            j += RESIDUE_LENGTH 
        i += RESIDUE_LENGTH
        
def center_coords(coords):
    av1 = sum(coords) / coords.shape[0]
    return coords - av1

def _superimpose(p1, p2, alignments, type):
    lmin = min(p1.shape[0], p2.shape[0])
    ##
    res = RESIDUE_LENGTH
    if lmin > 250:
        res = 50
    elif lmin > 200:
        res = 40
    elif lmin > 150:
        res = 30
    else:
        res = 20
    m1 = p1.shape[0] - res - 20
    m2 = p2.shape[0] - res - 20
    ## local structure alignment
    rs = []; rotrans = []; edgelist = []
    iter = 0
    for r1_idx in xrange(0, m1, res):
        for r2_idx in xrange(0, m2, res):
            iter += 1
            r1 = p1[r1_idx : r1_idx + res]
            r2 = p2[r2_idx : r2_idx + res]
            sup.set(r1, r2); sup.run(); u, t = sup.get_rotran()
            lp2 = np.dot(p2,u) + t
            edges = aligner[type](p1, lp2)
            sup.set(np.array(map(lambda x : p1[x[0]], edges)), np.array(map(lambda x : p2[x[1]], edges)))
            sup.run(); r = sup.get_rms(); u, t = sup.get_rotran()
            rs.append(r); rotrans.append((u, t)); edgelist.append(edges)
    print '%s' %[p1.shape[0], p2.shape[0], iter]
    ## find best overall alignment
    rs = np.array(rs)
    minindex = np.argmin(rs)
    rms = rs[minindex]
    rot, tran = rotrans[minindex]
    edges = edgelist[minindex]
    ## refine alignment with iterative DP
    lp2 = p2; pedges = edges
    for iter in xrange(5):
        sup.set(np.array(map(lambda x : p1[x[0]], edges)), np.array(map(lambda x : p2[x[1]], edges)))
        sup.run(); r = sup.get_rms(); u, t = sup.get_rotran()
        if r < rms:
            rms = r; rot = u; tran = t; pedges = edges
        lp2 = np.dot(p2, r) + t
        edges = aligner[type](p1, lp2) 
        
    return [rms, rot, tran, None, None, pedges]
  
def superimpose(p1, p2, alignments, type = 0):
    return _superimpose(p1, p2, alignments, type)
    
def ssuperimpose(p1, p2, alignments = None):
    return superimpose(p1, p2, alignments)

def nssuperimpose(p1, p2, alignments = None):
    return superimpose(p1, p2, alignments, 1)

##
data_dir = 'pdb_data'
file_names = os.listdir(data_dir) #[0:2]
#file_names = ['1aa9.pdb','1ash.pdb']
file_paths = map(lambda x : '%s/%s' %(data_dir, x), file_names)
pdb_parser = Bio.PDB.PDBParser(QUIET = True)
##
domain_structures = map(lambda x : pdb_parser.get_structure("d", x), file_paths)
ca_res_list = map(lambda x : get_ca_atom_list(x[0]), domain_structures)
domain_ca_atoms = map(lambda x : x[0], ca_res_list)
domain_res_atoms = map(lambda x : x[1], ca_res_list)
domain_coordinates = map(lambda x : np.array(map (lambda y: center_coords(y.get_coord()), x)), domain_ca_atoms)
for fn1, dc1 in zip(file_names, domain_coordinates):
    for fn2, dc2 in zip(file_names, domain_coordinates):
        if fn1 == fn2:
            continue
        [rms, rot, tran, p1_r, p2_r, edges] = ssuperimpose(dc1, dc2)
        print '%s %s %f %d' %(fn1, fn2, rms, len(edges))
        sys.stdout.flush()
    break

