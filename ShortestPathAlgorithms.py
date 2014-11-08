#!/usr/bin/python

"""
Implementations of two shortest path algorithms - bellman_ford, floyd_warshall.

Run this file as ./ShortestPathAlgorithms.py to get results for the required 
test cases. The test are coded into the main section of the source file
"""

import sys
INF = sys.maxint

class Vertex:
  def __init__(self, label):
    self.label = label
    self.distance = None
    self.predecessor = None
    
  def __repr__(self):
    return 'Label:%s, Distance:%s' %(self.label, self.distance)

class Edge:
  def __init__(self, source, destination, weight):
    self.source = source
    self.destination = destination
    self.weight = weight
    
  def __repr__(self):
    return 'Source:%s, Destination: %s, Weight: %d' %(self.source, self.destination, self.weight)
    
def bellman_ford(vertices, edges, source):
  for v in vertices:
    if v.label == source.label:
      v.distance = 0
    else:
      v.distance = INF
    v.predecessor = None
    
  for i in range(1, len(vertices) - 1):
    for uv in edges:
      u = uv.source
      v = uv.destination
      if u.distance + uv.weight < v.distance:
        v.distance = u.distance + uv.weight
        v.predecessor = u
        
  for uv in edges:
    u = uv.source
    v = uv.destination
    if u.distance + uv.weight < v.distance:
      return False
  
  return True    

def floyd_warshall(vertices, edges, source):
  ## init path
  path = {}
  for i in vertices:
    for j in vertices:
      if i == j:
        path[i,j] = 0
      else:
        path[i,j] = INF
  for uv in edges:
    path[uv.source,uv.destination] = uv.weight   
    
  for k in vertices:
    for i in vertices:
      for j in vertices:
        path[i,j] = min(path[i,j], path[i,k] + path[k,j])

  for k in vertices:
    if path[k,k] < 0:
      return False
  
  for v in vertices:
    v.distance = path[source,v]        
      
  return True

def do_test(function_name, vertices, edges, source):
  print '\tSource node label: %s' %source.label
  if function_name(vertices,edges,source) == True:
    print '\tDistances from nodes are:' 
    print '\t%s' %vertices
  else:
    print '\tGraph contains a negative-weight cycle'
      
if __name__ == '__main__':
  v1 = Vertex(1); v2 = Vertex(2); v3 = Vertex(3); v4 = Vertex(4)
  e1 = Edge(v1, v2, -3); e2 = Edge(v2, v3, 1); e3 = Edge(v3, v4, 2)
  e4 = Edge(v4, v1, 1); e5 = Edge(v2, v4, 2)
  vertices = [v1, v2, v3, v4]; edges = [e1, e2, e3, e4, e5]
  ## Test 1
  print 'Test case 1, Algorithm - Bellman-Ford'
  do_test(bellman_ford,vertices,edges,v2)
  ## Test 2
  print 'Test case 2, Algorithm - Bellman-Ford'
  e5.weight = 1
  do_test(bellman_ford,vertices,edges,v2)
  ## Test 3
  print 'Test case 3, Algorithm - Floyd-Warshall'
  e5.weight = 2
  do_test(floyd_warshall,vertices,edges,v2)
  ## Test 4
  print 'Test case 4, Algorithm - Floyd-Warshall'
  e5.weight = 1
  do_test(floyd_warshall,vertices,edges,v2)
  v1 = Vertex('r1');v2 = Vertex('r2');v3 = Vertex('r3');v4 = Vertex('r4');v5 = Vertex('r5');v6 = Vertex('r6')
  e1 = Edge(v1, v5, 2); e2 = Edge(v1, v4, 2); e3 = Edge(v1, v3, 3)
  e4 = Edge(v2, v1, 1); e5 = Edge(v2, v3, 1)
  e6 = Edge(v3, v4, -1)
  e7 = Edge(v5, v3, -1); e8 = Edge(v5, v4, -2)
  e9 = Edge(v6,v1,0);e10 = Edge(v6,v2,0);e11 = Edge(v6,v3,0);e12 = Edge(v6,v4,0);e13 = Edge(v6,v5,0);
  vertices = [v1, v2, v3, v4, v5, v6]
  edges = [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13]
  ## Test 5
  print 'Test case 5, Algorithm - Bellman-Ford'
  do_test(bellman_ford,vertices,edges,v6)
  ## Test 6
  print 'Test case 6, Algorithm - Floyd-Warshall'
  do_test(floyd_warshall,vertices,edges,v6)
  ## Test 7
  print 'Test case 7, Algorithm Bellman-Ford'
  do_test(bellman_ford,vertices,[],v2)
  ## Test 7
  print 'Test case 8, Algorithm Floyd-Warshall'
  do_test(floyd_warshall,vertices,[],v2)

