#!/usr/bin/python

import Bio.PDB
import numpy as np
import sys
from datetime import datetime
from numpy import sum, sqrt, dot
from Bio.PDB.QCPSuperimposer import QCPSuperimposer
from Bio.PDB.PSC.align import Align
from scipy.spatial import distance
import multiprocessing

RESIDUE_LENGTH = 20
RMS_THRESH = 0.001
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

def _task(task):
    p1, p2, i, j, type = task
    sup.set(p1[i:i+RESIDUE_LENGTH], p2[j:j+RESIDUE_LENGTH]); sup.run()  
    u, t = sup.get_rotran()
    lp2 = dot(p2, u) + t
    edges = aligner[type](p1, lp2)
    a = np.array(map(lambda x : p1[x[0]], edges))
    b = np.array(map(lambda x : lp2[x[1]], edges))
    sup.set(a, b); sup.run()
    r = sup.get_rms()
    return [i,j,(u,t),edges,r]
  
def _superimpose(p1, p2, type):
    tasks = []; i = 0
    while i < p1.shape[0] - RESIDUE_LENGTH:
        j = 0
        while j < p2.shape[0] - RESIDUE_LENGTH:
            tasks.append([p1,p2,i,j,type])
            j += RESIDUE_LENGTH
        i += RESIDUE_LENGTH
    for number_processes in range(1, 16):
        pool = multiprocessing.Pool(number_processes)
        t0 = datetime.now()
        UsAndTs = pool.map(_task, tasks)
        dif = datetime.now()-t0
        if number_processes == 1:
            told = dif.total_seconds() * 1000
            tnew = told
        else:
            tnew = dif.total_seconds() * 1000
        speedup = told * 1.0 / tnew
        efficiency = speedup / number_processes
        print '%d %f %f' %(number_processes,speedup,efficiency)
        
    ###
    ret = np.array(map(lambda x: x[4], UsAndTs))
    i, j, (rot, tran), pedges, rms = UsAndTs[np.argmin(ret)]
    return [rms, rot, tran, None, None, pedges]
    
def superimpose(p1, p2, type = 0):
    t0 = datetime.now()
    [rms, rot, tran, p1_r, p2_r, edges] = _superimpose(ref_coords, sample_coords, type)
    dif = datetime.now()-t0
    print "#%s: %d, %s, %d (msec)" %(type, len(edges), [len(ref_atoms), len(sample_atoms)], dif.total_seconds() * 1000)
    print "#%s" %([('%s_%d' %(ref_reses[i].get_resname(), i), '%s-%d' %(sample_reses[j].get_resname(),j)) for i,j in edges])
    print "#rmsd: %f" %rms
    print "#rot: %s" %(','.join(map(lambda x: '%f %f %f' %(x[0],x[1],x[2]), rot)))
    print "#tran: %s" %tran
    
def ssuperimpose(p1, p2):
    superimpose(p1, p2)

def nssuperimpose(p1, p2):
    superimpose(p1, p2, 1)

## PDB domains
if len(sys.argv) < 3:
   print 'usage: python test_psalign.py dom1 dom2 '
   sys.exit(-1)
pdomain1 = '%s.pdb' %sys.argv[1]
pdomain2 = '%s.pdb' %sys.argv[2]
print '#%s %s' %(pdomain1,pdomain2)
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
#nssuperimpose(ref_coords, sample_coords)                                                                                                    
