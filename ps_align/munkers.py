#!/usr/bin/python

import munkers
import numpy as np

matrix = np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]]).astype(np.int32)
print matrix
print np.array(munkers.run_munkers(matrix, 1)).reshape(matrix.shape)

print "done"