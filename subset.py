#!/usr/bin/python
#
# Follow naming conventions described at:
#  PEP 0008 (https://www.python.org/dev/peps/pep-0008/)
#

import numpy as np

class SubsetGenerator:
    def __init__(self, fname=None):
        self._fname = fname
        self._FVMIN = 0
        self._SEP1 = ':'
        self._SEP2 = ' '
        self._skipf = lambda x : not x.startswith('#')
        if fname != None:
            self.data = self.__readsparsedata()

    def __read_line(self, l):
        zipped = map(lambda x : tuple(x.split(self._SEP1)), l.split(self._SEP2))
        idx, data = zip(*zipped)
        idx = map(int, idx); data = map(float, data)
        ret = np.zeros(max(idx) + 1)
        ret[idx] = data
        return ret

    def __read_sparse_data(self):
        # read data
        vectors = map(self.__readline, filter(self._skipf, open(self._fname)))
        # find longest vector
        max_entries = max(map(len, vectors))
        # create masked vectors
        maskf = lambda b : np.ma.array(np.resize(b,max_entries),
            mask=np.concatenate([np.zeros(len(b),dtype=bool),
                                 np.ones(max_entries-len(b), dtype=bool)]))
        # return square np matrix
        return np.array(map(maskf, data))

    def get_data(self): return self.data
        
        
if __name__ == '__main__':
    test_file_name = 'subset.dat'
    s = SubsetGenerator(test_file_name)
    print s.get_data()
