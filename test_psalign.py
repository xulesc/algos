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
  
def ssuperimpose(p1, p2):
    sup = QCPSuperimposer(); align = Align()
    l1 = p1.shape[0]; l2 = p2.shape[0]; scores = []
    lp2 = p2; min_rot = []; min_tran = []; old_rms = sys.maxint
    for iter in xrange(10):
        min_rms = sys.maxint; 
        for x in range(0, l1 - RESIDUE_LENGTH):
            r1 = p1[x : x + RESIDUE_LENGTH]; 
            for y in range(0, l2 - RESIDUE_LENGTH):
                r2 = lp2[y : y + RESIDUE_LENGTH]; 
                sup.set(r1, r2)
                sup.run()
                lrms = sup.get_rms()
                if lrms < min_rms:
                    min_rms = lrms
                    min_rot, min_tran = sup.get_rotran()
# optimal from tmalign for 1aa9 & 1ash
#        min_tran = np.array([-27.1331353550,-147.9265112953,0.9599895782])              
#        min_rot  = np.array([
#            [0.3141317101,0.2302470691,-0.9210361317],
#            [0.9493282736,-0.0862496500,0.3022198320],
#            [-0.0098538134,-0.9693024735,-0.2456738027]
#        ])
        pr2 = dot(lp2, min_rot) + min_tran
        edges = align.sequential(p1, pr2)
        a = np.array(map(lambda x : p1[x[0]], edges));  b = np.array(map(lambda x : pr2[x[1]], edges));
        d = a - b
        arms = sqrt(dot(d,d.T).diagonal().sum()/len(edges))
        if arms >= old_rms:
            return [old_rms, old_rot, old_tran, [], [], []]
        lp2 = pr2
        old_rms = arms
        old_rot = min_rot
        old_tran = min_tran

def nssuperimpose(p1, p2):
    frms = 0; frot = []; ftran = []; p1_r = []; p2_r = []; medges = 0
    mrms = sys.float_info.max
    l1 = p1.shape[0]; l2 = p2.shape[0]; ml = max(l1, l2)
    sup = QCPSuperimposer()
    align = Align()
    for x in range(0, l1 - RESIDUE_LENGTH, RESIDUE_LENGTH):
        r1 = p1[x : x + RESIDUE_LENGTH]
        for y in range(0, l2 - RESIDUE_LENGTH, RESIDUE_LENGTH):
            r2 = p2[y : y + RESIDUE_LENGTH]
            ##
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
            pr2 = dot(p2, rot) + tran
            edges = align.non_sequential(p1, pr2)
            a = np.array(map(lambda x : p1[x[0]], edges));  b = np.array(map(lambda x : pr2[x[1]], edges));  
            d = a - b
            arms = sqrt(dot(d,d.T).diagonal().sum()/len(edges))
            if arms < RMS_THRESH:
                return [arms, rot, tran, [x, x + RESIDUE_LENGTH-1], [y, y + RESIDUE_LENGTH-1], edges]
            if mrms > arms:
                mrms = arms
                frot = rot
                ftran = tran
                p1_r = [x, x + RESIDUE_LENGTH-1]
                p2_r = [y, y + RESIDUE_LENGTH-1]
                medges = edges
    return [mrms, frot, ftran, p1_r, p2_r, medges]
            
## PDB domains
pdomain1 = '1aa9.pdb'; pdomain2 = '1ash.pdb'
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
##
## non-sequential
t0 = datetime.now()
[rms, rot, tran, p1_r, p2_r, edges] = nssuperimpose(ref_coords, sample_coords)
dif = datetime.now()-t0
print "non-sequential threading: %d, %s, %d (msec)" %(len(edges), [len(ref_atoms), len(sample_atoms)], dif.total_seconds() * 1000)
#print [('%s_%d' %(ref_reses[i].get_resname(), i), '%s-%d' %(sample_reses[j].get_resname(),j)) for i,j in edges]
print "rmsd: %f" %rms
print "rot: %s" %rot
print "tran: %s" %tran
#print "residues: %s:%s"  %(p1_r, p2_r)
## sequential
t0 = datetime.now()
[rms, rot, tran, p1_r, p2_r, edges] = ssuperimpose(ref_coords, sample_coords)
dif = datetime.now()-t0
print "sequential alignemnt: %d (msec)" %(dif.total_seconds() * 1000)
print "rmsd: %f" %rms
print "rot: %s" %rot
print "tran: %s" %tran
#print "residues: %s:%s"  %(p1_r, p2_r)
                                                                                                    
