#!/usr/bin/env python

import numpy as np
from scipy.spatial import distance
import Bio.PDB
import zlib
import time
import sys
from collections import defaultdict
from collections import Counter
import itertools
from scipy.spatial.distance import cdist

RESIDUE_LENGTH=8

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
	
def make_residues2(coords): return map(lambda x : '%s' %center_coords(coords[x:x+RESIDUE_LENGTH]).tolist(), range(0,len(coords)))

def make_rep(res):
	return sorted(set(cdist(res,res,'euclidean').flatten()))

def make_residues(coords):
	return map(lambda x : '%s' %make_rep(center_coords(coords[x:x+RESIDUE_LENGTH])), range(0,len(coords)))

def center_coords(coords): return coords - (sum(coords) / len(coords))
	
def make_ngram_hash(inverted_index, dom_coords, dom_coords_idx):
	coords = center_coords(dom_coords)
	return map(lambda x : inverted_index[x].append(dom_coords_idx), make_residues(coords))
	
def do_query(dom_coords,inverted_index):
	residues = make_residues(center_coords(dom_coords))
	l2 = map(lambda x : inverted_index[x], residues)
	dom_indices = list(itertools.chain(*l2))
	print Counter(dom_indices)

data_dir = '../pdb_data'; file_names = ['6xia.pdb','1ash.pdb','1aa9.pdb','1bab.pdb','1cd8.pdb']
file_paths = map(lambda x : '%s/%s' %(data_dir, x), file_names)
pdb_parser = Bio.PDB.PDBParser(QUIET = True)
## load domains and extract coordinates
domain_structures = map(lambda x : pdb_parser.get_structure("d", x), file_paths)
ca_res_list = map(lambda x : get_ca_atom_list(x[0]), domain_structures)
domain_ca_atoms = map(lambda x : x[0], ca_res_list)
domain_coordinates = map(lambda x : np.array(map (lambda y: y.get_coord(), x)), domain_ca_atoms)
## make inverted index
inverted_index = defaultdict(list)
for i in range(0,len(file_names)):
	make_ngram_hash(inverted_index,domain_coordinates[i],i)
##print inverted_index
## query inverted index
for i in range(0,len(file_names)):
	do_query(domain_coordinates[i],inverted_index)

