#!/usr/bin/python

import munkers
import numpy as np

matrix = np.arange(25).reshape(5,5).astype(np.int32)
munkers.run_munkers(matrix, 1)

print "done"