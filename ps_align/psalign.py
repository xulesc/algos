#!/usr/bin/python

import Bio.PDB
import munkers
import numpy as np
from scipy.spatial import distance
import sys
from datetime import datetime
from numpy import sum, sqrt, dot
from Bio.QCPSuperimposer import QCPSuperimposer

RESIDUE_LENGTH = 8
RMS_THRESH = 0.00001

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
  
def superimpose(p1, p2):
    frms = 0; frot = []; ftran = []; p1_r = []; p2_r = []
    mrms = sys.float_info.max
    l1 = p1.shape[0]; l2 = p2.shape[0]
    for x in range(0, l1 - RESIDUE_LENGTH, RESIDUE_LENGTH):
        r1 = p1[x : x + RESIDUE_LENGTH]
        for y in range(0, l2 - RESIDUE_LENGTH, RESIDUE_LENGTH):
            r2 = p2[y : y + RESIDUE_LENGTH]
            ##
            sup = QCPSuperimposer()
            sup.set(r1,r2)
            sup.run()
            rms = sup.get_rms()
            rot, tran = sup.get_rotran()
            ###################################################################################################################
            ## Align using bipartite graph techniques based on description given in Wang, Y.; Makedon, F.; Ford, J.; Heng Huang, 
            ## "A bipartite graph matching framework for finding correspondences between structural elements in two proteins," 
            ## Engineering in Medicine and Biology Society, 2004. IEMBS '04. 26th Annual International Conference of the IEEE , 
            ## vol.2, no., pp.2972,2975, 1-5 Sept. 2004
            ## Note we are not rotating anything just finding a global non-sequential alignment
            ###################################################################################################################
            dist_matrix = distance.cdist(p1, dot(p2, rot) + tran, 'euclidean').astype(np.int32)
            cost_matrix = np.array(munkers.run_munkers(dist_matrix, 0)).reshape(dist_matrix.shape)
            non_zero = cost_matrix > 0
            edges = np.column_stack(np.where(non_zero))
            a = np.array(map(lambda x : p1[x[0]], edges));  b = np.array(map(lambda x : p2[x[1]], edges));  
            d = a - b
            arms = sqrt(sum(sum(d*d))/len(edges))
            sys.stdout.flush()
            if arms < RMS_THRESH:
                return [arms, rot, tran]
            if mrms > arms:
                mrms = arms
                frot = rot
                ftran = tran
                p1_r = [x, x + RESIDUE_LENGTH-1]
                p2_r = [y, y + RESIDUE_LENGTH-1]
    return [frms, frot, ftran, p1_r, p2_r]
            
## PDB domains
pdomain1 = '1aa9.pdb'; pdomain2 = '6xia.pdb'
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
##
[rms, rot, tran, p1_r, p2_r] = superimpose(ref_coords, sample_coords)
dist_matrix = distance.cdist(ref_coords, dot(sample_coords, rot) + tran, 'euclidean').astype(np.int32)
t0 = datetime.now()
cost_matrix = np.array(munkers.run_munkers(dist_matrix, 0)).reshape(dist_matrix.shape)
non_zero = cost_matrix > 0
edges = np.column_stack(np.where(non_zero))
dif = datetime.now() - t0
print "non-sequential threading: %d, %s, %d (msec)" %(len(edges), dist_matrix.shape, dif.total_seconds() * 1000)
print [('%s_%d' %(ref_reses[i].get_resname(), i), '%s-%d' %(sample_reses[j].get_resname(),j)) for i,j in edges]
print "rmsd: %f" %rms
print "rot: %s" %rot
print "tran: %s" %tran
print "residues: %s:%s"  %(p1_r, p2_r)
                                                                                                    
