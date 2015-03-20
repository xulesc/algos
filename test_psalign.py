#!/usr/bin/python

import Bio.PDB
import numpy as np
import sys, os
from datetime import datetime
from numpy import sum, sqrt, dot
from Bio.PDB.QCPSuperimposer import QCPSuperimposer
from Bio.PDB.PSC.align import Align
from scipy.spatial import distance
import multiprocessing
import random

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

#------------ Rotation matrix to rotate Chain_1 to Chain_2 ----------
# m          t(m)         u(m,1)         u(m,2)         u(m,3)
# 1    -27.1331353550   0.3141317101   0.2302470691  -0.9210361317
# 2   -147.9265112953   0.9493282736  -0.0862496500   0.3022198320
# 3      0.9599895782  -0.0098538134  -0.9693024735  -0.2456738027
#        u = np.array([[0.3141317101,0.2302470691,-0.9210361317],
#          [0.9493282736,-0.0862496500,0.3022198320],
#          [-0.0098538134,-0.9693024735,-0.2456738027]])
#        t = np.array([-27.1331353550,-147.9265112953,0.9599895782])
def _superimpose(p1, p2, type):
    lmin = min(p1.shape[0], p2.shape[0])
    rms = sys.maxint; prms = -1;
    while True:
        if abs(rms-prms) < RMS_THRESH:
             break
        prms = rms; lp2 = p2
        #
        i1 = random.sample(range(0, lmin), lmin/2)
        i2 = random.sample(range(0, lmin), lmin/2)
        i1.sort(); i2.sort()
        edges = zip(i1, i2)
        #
        for loop in xrange(5):        
             a = np.array(map(lambda x : p1[x[0]], edges))
             b = np.array(map(lambda x : lp2[x[1]], edges))
             sup.set(a, b); sup.run()
             r = sup.get_rms() / len(edges)
             u, t = sup.get_rotran()
             if r < rms:
                 rms = r
                 rot = u
                 tran = t
                 pedges = edges
             lp2 = np.dot(lp2, u) + t                        
             edges = aligner[type](p1, lp2)         
    return [rms, rot, tran, None, None, pedges]
  
def _psuperimpose(p1, p2, type):
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
#    t0 = datetime.now()
    #[rms, rot, tran, p1_r, p2_r, edges] = 
    return _superimpose(p1, p2, type)
#    dif = datetime.now()-t0
#    print "#%s: %d, %s, %d (msec)" %(type, len(edges), [len(ref_atoms), len(sample_atoms)], dif.total_seconds() * 1000)
#    print "#%s" %([('%s_%d' %(ref_reses[i].get_resname(), i), '%s-%d' %(sample_reses[j].get_resname(),j)) for i,j in edges])
#    print "#rmsd: %f" %rms
#    print "#rot: %s" %(','.join(map(lambda x: '%f %f %f' %(x[0],x[1],x[2]), rot)))
#    print "#tran: %s" %tran
    
    
def ssuperimpose(p1, p2):
    return superimpose(p1, p2)

def nssuperimpose(p1, p2):
    return superimpose(p1, p2, 1)

##
data_dir = 'pdb_data'
file_names = os.listdir(data_dir)[0:2]
file_paths = map(lambda x : '%s/%s' %(data_dir, x), file_names)
pdb_parser = Bio.PDB.PDBParser(QUIET = True)
domain_structures = map(lambda x : pdb_parser.get_structure("d", x), file_paths)
domain_ca_atoms = map(lambda x : get_ca_atom_list(x[0])[0], domain_structures)
domain_coordinates = map(lambda x : np.array(map (lambda y: y.get_coord(), x)), domain_ca_atoms)
for fn1, dc1 in zip(file_names, domain_coordinates):
    for fn2, dc2 in zip(file_names, domain_coordinates):
        [rms, rot, tran, p1_r, p2_r, edges] = ssuperimpose(dc1, dc2)
        print '%s %s %f' %(fn1, fn2, rms)
        sys.stdout.flush()

