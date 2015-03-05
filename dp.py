#!/usr/bin/python

import numpy as np

def dp(gapcost, S):
    w = gapcost; (rows, cols) = S.shape
    M = np.zeros((rows + 1, cols + 1), dtype=np.int32)
    ## init
    for i in xrange(rows + 1):
        M[i][0] = 0
    for j in xrange(cols + 1):
        M[0][j] = 0
    ## fill 
    for i in xrange(1, rows + 1):
        for j in xrange(1, cols + 1):
            M[i][j] = max(M[i-1][j-1]+S[i-1][j-1], M[i][j-1]+w, M[i-1][j]+w)
    print M
    #idxs = np.where(M == M.max())
    #(X, Y) = np.where(M == M.max()); X = X[len(X)-1]; Y = Y[len(Y)-1]
    #print '%d %d' %(X,Y)
    X = rows - 1; Y = cols - 1
    ## traceback
    path = []
    path.append((X, Y))
    while X >= 0 and Y >= 0:
        current = M[X][Y]
        diag = M[X-1][Y-1]
        up = M[X-1][Y]
        left = M[X][Y-1]
        if diag >= up and diag >= left:
            print '%s' %[diag,up,left]
            path.append((X, Y))
            X -= 1
            Y -= 1
        elif up >= diag and up >= left:
            X -= 1
        else:
            Y -= 1
    path.reverse()
    return path    

if __name__ == '__main__':
    rows = 7; cols = 11; gapcost = -2
    # S 1 for match, 0 for mismatch
    S = [[1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 
    [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 
    [0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1], 
    [0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0], 
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], 
    [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 
    [0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1]]    
    
    S = [[2, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1], 
    [2, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1], 
    [-1, 2, 2, -1, -1, -1, 2, -1, -1, -1, 2], 
    [-1, -1, -1, 2, 2, -1, -1, -1, 2, 2, -1], 
    [-1, -1, -1, -1, -1, 2, -1, -1, -1, -1, -1], 
    [2, -1, -1, -1, -1, -1, -1, 2, -1, -1, -1], 
    [-1, 2, 2, -1, -1, -1, 2, -1, -1, -1, 2]]    
    
    S = np.array(S)    
    print dp(gapcost, S)
