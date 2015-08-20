#!/usr/bin/python
#
# Follow naming conventions described at:
#  PEP 0008 (https://www.python.org/dev/peps/pep-0008/)
#

import numpy as np
import random

class SubsetGenerator:

    """ Return a representative subset of a given matrix.

    Generate row based subset from a given data matrix.

    The class provides functionality for selecting rows from a sparse matrix of
    feature vectors. The subset generated is expected to be a representative of
    the original matrix making use of the sampling method described in 
    http://www.cs.rpi.edu/~drinep/Papers/Drineas_MFO_04.pdf.

    Input data to the class is a sparse matrix of the data however internally
    it is stored as a dense matrix also the subset generated is returned as a 
    dense Numpy matrix.
    """

    def __init__(self, fname=None):
        self._fname = fname
        self._FVMIN = 0
        self._SEP1 = ':'
        self._SEP2 = ' '
        self._skipf = lambda x : not x.startswith('#')
        if fname != None:
            self.data = self.__read_sparse_data()

    def __read_line(self, l):
        zipped = map(lambda x : tuple(x.split(self._SEP1)), l.split(self._SEP2))
        idx, data = zip(*zipped)
        idx = map(int, idx); data = map(float, data)
        ret = np.zeros(max(idx) + 1)
        ret[idx] = data
        return ret

    def __read_sparse_data(self):

        """ Load spare matrix data to dense Numpy array.

        The method loads sparse data from a file into a dense array. The data
        is expected to be in the following format e.g.

        0:.56 10:5.6
        1:8.2 2:3.2
        ...

        Matrix rows are expected to be separated by a newline and the columns
        are separated by a space. Each column data is expected to be in the 
        format <column_index>:<column_value>. Rows starting with # are excluded
        as comment lines.

        Minimum column_value is assumed to be 0 and missing values are filled
        in with 0s.
        """
        
        # read data
        vectors = map(self.__read_line, filter(self._skipf, open(self._fname)))
        # find longest vector
        max_entries = max(map(len, vectors))
        # create masked vectors
        maskf = lambda b : np.ma.array(np.resize(b,max_entries),
            mask=np.concatenate([np.zeros(len(b),dtype=bool),
                                 np.ones(max_entries-len(b), dtype=bool)]))
        # return dense np matrix
        return np.array(map(maskf, vectors))

    def get_data(self): return self.data

    def set_data(self, data): self.data = data

    def make_subset(self, similarity=0.25):

        """ Return a subset of data matrix.

        The method generates a subset of the data matrix as described in the
        aforementioned paper.

        L2-norms are calculated for all rows of the original data matrix. The
        norms are converted to probabilities by dividing with the sum of the
        norm vector. Uniform sampling from the cumulative sum vector over the
        norm-probabilitiy vector is used for generating the subset. Stopping
        criteria for the subset generation is the similarity between the sum
        of norms of all rows of the original data matrix and the subset matrix.

        The similarity argument accepted by the method has a default value of
        0.25 i.e. 25% similarity between the sum of norms of the subset and 
        the original data matrix.
        """
        
        # make L2-norm per row
        norms = map(np.linalg.norm, self.data)
        # convert norms to probabilities
        norms /= sum(norms)
        # make CDF
        cdf = np.cumsum(norms)
        # make subset
        visited = {}
        # loop till similarity criteria is met. using map handles duplicate
        # row selection.
        while sum(norms[visited.keys()]) < similarity:
            visited[np.searchsorted(cdf, random.uniform(0, 1), sorter=None)] = 1
        return self.data[visited.keys()]
        
if __name__ == '__main__':
    test_file_name = 'subset.dat'
    s = SubsetGenerator(test_file_name)
    print s.get_data()
    print s.make_subset()
