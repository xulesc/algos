#!/usr/bin/python

import numpy as np
from scipy.spatial import distance
import Bio.PDB
import zlib
import time
import sys

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

def lzw(instring = None):
	return [None, zlib.compress(instring)]

def _make_contact_map(coords):
	dist = distance.cdist(coords, coords, 'euclidean')
	idx1, idx2 = np.where(dist <= 6.5)
	fidx = np.where((idx2 - idx1) >= 2)
	return [idx1[fidx], idx2[fidx], len(fidx[0]), len(coords)]

def cmap2str(contacts = None):
	if contacts == None:
		return ''
	return '.'.join(map(lambda x : '%s.%s' %(x[0],x[1]), contacts))
			
def get_contact_map(coords):
	[X, Y, c, a] = _make_contact_map(coords[0:min(len(coords), 1000)])
	return zip(X, Y)

def compress(contacts):
	t0 = time.clock()
	s = cmap2str(contacts)
	print time.clock() - t0, " seconds tostr time"
	t0 = time.clock()
	[d1, ot1] = lzw(s)
	print time.clock() - t0, " seconds compress time"
	return ot1

if __name__ == '__main__':
	#print _compress('bannana_bandana')
	#import sys;sys.exit(-1);
	data_dir = '../pdb_data'; file_names = ['6xia.pdb','6xia.pdb']
	file_paths = map(lambda x : '%s/%s' %(data_dir, x), file_names)
	pdb_parser = Bio.PDB.PDBParser(QUIET = True)
	## load domains and extract coordinates
	domain_structures = map(lambda x : pdb_parser.get_structure("d", x), 
		file_paths)
	ca_res_list = map(lambda x : get_ca_atom_list(x[0]), domain_structures)
	domain_ca_atoms = map(lambda x : x[0], ca_res_list)
	domain_coordinates = map(lambda x : np.array(
		map (lambda y: y.get_coord(), x)), domain_ca_atoms)
	## get contact maps
	t0 = time.clock()
	cm1 = get_contact_map(domain_coordinates[0])
	print time.clock() - t0, " seconds cmap time"
	t0 = time.clock()
	cm2 = get_contact_map(domain_coordinates[1])
	print time.clock() - t0, " seconds cmap time"
	##
	ot1 = compress(cm1); ot2 = compress(cm2)
	ot3 = compress(cm1+cm2); ot4 = compress(cm2+cm1)
	x = 1.0 * len(ot1); y = 1.0 * len(ot2); 
	xy = 1.0 * len(ot3); yx = 1.0 * len(ot4)
	print 'dist: %f' %(max(yx-y,xy-x)/max(y,x))

