import Bio, Bio.PDB
import numpy as np
import time
from scipy.spatial import distance
import matplotlib.pylab as plt
from numpy import linalg as LA

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

def make_contact_map(coords):
	dist = distance.cdist(coords, coords, 'euclidean')
	ret = np.zeros(dist.shape)
	## "protein structure comparison vis contact map alignment" 
	## (felipe leal valentim)
	ret[(dist >= 5) & (dist <= 12)] = 1
	return ret

def show_contact_maps(cm1,cm2):
	fig = plt.figure()
	ax = fig.add_subplot(2,1,1)
	ax.set_aspect('equal')
	plt.imshow(cm1, interpolation='nearest', cmap=plt.cm.ocean)
	plt.colorbar()
	ax = fig.add_subplot(2,1,2)
	ax.set_aspect('equal')
	plt.imshow(cm2, interpolation='nearest', cmap=plt.cm.ocean)
	plt.colorbar()
	plt.show()	

if __name__ == '__main__':
	data_dir = '../pdb_data'; file_names = ['1aa9.pdb','1ash.pdb']
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
	cm1 = make_contact_map(domain_coordinates[0])
	print time.clock() - t0, " seconds cmap time"
	t0 = time.clock()
	cm2 = make_contact_map(domain_coordinates[1])
	print time.clock() - t0, " seconds cmap time"
	print cm1.diagonal().sum()
	print cm2.diagonal().sum()
	## visualize the contact maps
	if True:
		show_contact_maps(cm1, cm2)
	## eigen decomposition
	w1, v1 = LA.eig(cm1)	
	w2, v2 = LA.eig(cm2)
	print [w1.shape, w2.shape]

