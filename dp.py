#!/usr/bin/python

import numpy as np

def dp2(gapcost, S):
    w = gapcost; (rows, cols) = S.shape
    ## init
    dir = np.zeros((rows + 1, cols + 1), dtype=np.int32)
    val = np.zeros((rows + 1, cols + 1), dtype=np.int32)
    ## decide matrix and path
    for i in xrange(1, rows):
        for j in xrange(1, cols):
            D = val[i-1][j-1] + S[i-1][j-1]
            H = val[i-1][j]
            if dir[i-1][j] == 1:
                H = H + gapcost
            V = val[i][j-1]
            if dir[i][j-1] == 1:
                V = V + gapcost
            
            if D >= H and D >= V:
                dir[i][j] = 1
                val[i][j] = D
            else:
                dir[i][j] = 0
                val[i][j] = max(V, H)
    print dir
    print val
    ## extract the alignment
    (i, j) = (rows, cols); path = [(i - 1, j - 1)]
    while i > 0 and j > 0:
        if dir[i][j] == 1:
            path.append((i - 1, j - 1))
            i -= 1
            j -= 1
        else:
            H = val[i-1][j]
            if dir[i-1][j] == 1:
                H = H + gapcost
            V = val[i][j-1]
            if dir[i][j-1] == 1:
                V = V + gapcost
            if V >= H:
                j -= 1
            else:
                i -= 1
    path.reverse()
    return path
  
  
#http://www.avatar.se/molbioinfo2001/dynprog/adv_dynamic.html
def dp(gapcost, S):
    w = gapcost; (rows, cols) = S.shape
    # init
    M = np.zeros((rows + 1, cols + 1), dtype=np.float)
    ## fill 
    for i in xrange(1, rows + 1):
	for j in xrange(1, cols + 1):
	    M[i][j] = max(M[i-1][j-1]+S[i-1][j-1], M[i][j-1]+w, M[i-1][j]+w)
    X = rows; Y = cols
    ## traceback
    alignment = []
    while X > 0 and Y > 0:
        ## diagonal, left, up
	u = [(X-1,Y-1),(X-1,Y),(X,Y-1)]
	v = np.array([M[X-1][Y-1] + S[X-1][Y-1], M[X-1][Y] + w, M[X][Y-1] + w])
	## which value is equal to the current value
	i = np.where(v==M[X][Y])[0][0]
	(X,Y) = u[i]
	if i == 0:
	    alignment.append((X, Y))
    alignment.reverse()
    return alignment    

if __name__ == '__main__':
    rows = 7; cols = 11; 
    # S 1 for match, 0 for mismatch
    gapcost = 0
    S = [[1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 
    [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 
    [0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1], 
    [0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0], 
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], 
    [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 
    [0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1]]    
    
    gapcost = -2
    S = [[2, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1], 
    [2, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1], 
    [-1, 2, 2, -1, -1, -1, 2, -1, -1, -1, 2], 
    [-1, -1, -1, 2, 2, -1, -1, -1, 2, 2, -1], 
    [-1, -1, -1, -1, -1, 2, -1, -1, -1, -1, -1], 
    [2, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1], 
    [-1, 2, 2, -1, -1, -1, 2, -1, -1, -1, 2]]    
    
    S = np.array(S)    
    print dp(gapcost, S)
