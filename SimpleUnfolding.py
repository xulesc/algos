#!/usr/bin/python

"""
Unfoling dsp implementations
"""

import math

class Node:
  def __init__(self, name):
    self.name = name
  def __repr__(self):
    return '%s' %self.name
    
class Edge:
  def __init__(self, source, target, delay = 0):
    self.source = source
    self.target = target
    self.delay = delay
  def __repr__(self):
    return 'Source: %s: Target: %s: Delay: %sD' %(self.source, self.target, self.delay)
    
def simple_unfold(nodes, edges, unfold):
  new_edges = []
  for edge in edges:
    for i in range(0, unfold):
      source = Node('%s%d' %(edge.source.name, i))      
      target = Node('%s%d' %(edge.target.name, ((i + edge.delay) % unfold)))
      delay = int(math.floor((i + edge.delay) / unfold))
      new_edge =  Edge(source, target, delay)
      new_edges.append(new_edge)
  return new_edges

def run(nodes, edges, unfold, graph_number, print_output = True):
  new_edges = simple_unfold(nodes, edges, unfold)
  if print_output == True:
    console_display(new_edges)
  return graph_as_string(new_edges, graph_number)

def console_display(edges):
  for edge in edges:
    print '\t%s' %edge

def graph_as_string(edges, number):
  ret = 'digraph g%d {\nnode [shape=circle]' %number
  for edge in edges:
    source = edge.source.name
    target = edge.target.name
    delay = edge.delay
    if delay > 0:
      label = '\"%dD\"' %delay
      ret = ret + '%s -> %s [label=%s]\n' %(source,target,label)      
    else:
      ret = ret + '%s -> %s \n' %(source,target)
  ret = ret + '}\n'
  return ret
  
if __name__ == '__main__':
  graph_string = ''
  ## Ex 1
  na = Node('A'); nb = Node('B'); nc = Node('C'); nd = Node('D')
  e1 = Edge(na, nc); e2 = Edge(nc, nd, 9); e3 = Edge(nc, nb); e4 = Edge(nd, nc)
  nodes = [na, nb, nc, nd]; edges = [e1, e2, e3, e4]
  graph_string = graph_string + run(nodes, edges, 2, 1)
  
  ## Ex 2
  na = Node('A'); nb = Node('B'); nc = Node('C')
  e1 = Edge(na,nb,10); e2 = Edge(nb,nc,2); e3 = Edge(nc,na,20)
  nodes = [na,nb,nc]; edges = [e1,e2,e3]
  graph_string = graph_string + run(nodes,edges,3, 2)
  #  
  na = Node('A'); nb = Node('B'); nc = Node('C')
  e1 = Edge(na,nb,10); e2 = Edge(nb,nc,2); e3 = Edge(nc,na,20)
  nodes = [na,nb,nc]; edges = [e1,e2,e3]
  graph_string = graph_string + run(nodes,edges,4, 3)
    
  ## Ex 3
  na = Node('A'); nb = Node('B'); nc = Node('C')
  e1 = Edge(na,nb,10); e2 = Edge(nb,nc,3); e3 = Edge(nc,na,16)
  nodes = [na,nb,nc]; edges = [e1,e2,e3]
  graph_string = graph_string + run(nodes,edges,2,4)
  #
  na = Node('A'); nb = Node('B'); nc = Node('C')
  e1 = Edge(na,nb,10); e2 = Edge(nb,nc,3); e3 = Edge(nc,na,16)
  nodes = [na,nb,nc]; edges = [e1,e2,e3]
  graph_string = graph_string + run(nodes,edges,5,5)
  #
  na = Node('A'); nb = Node('B'); nc = Node('C'); nd = Node('D'); ne = Node('E')
  e1 = Edge(na,nb,1); e2 = Edge(nb,nc); e3 = Edge(nc,ne,7); e4 = Edge(ne, nb, 3); e5 = Edge(ne,nd,2); e6=Edge(nd,na)
  nodes = [na,nb,nc,nd,ne]; edges = [e1,e2,e3,e4,e5,e6]
  graph_string = graph_string + run(nodes,edges,2,6)
  #
  na = Node('A'); nb = Node('B'); nc = Node('C'); nd = Node('D'); ne = Node('E')
  e1 = Edge(na,nb,1); e2 = Edge(nb,nc); e3 = Edge(nc,ne,7); e4 = Edge(ne, nb, 3); e5 = Edge(ne,nd,2); e6=Edge(nd,na)
  nodes = [na,nb,nc,nd,ne]; edges = [e1,e2,e3,e4,e5,e6]
  graph_string = graph_string + run(nodes,edges,5,7)
    
  ## Ex 4
  na = Node('A'); nb = Node('B'); nc = Node('C'); nd = Node('D'); ne = Node('E'); nf = Node('F'); ng = Node('G')
  e1=Edge(na,nb,2);e2=Edge(nb,nc,3);e3=Edge(na,nf);e4=Edge(nb,nd);e5=Edge(nd,nf);e6=Edge(nc,ne);
  e7=Edge(nf,ng);e8=Edge(ne,ng)
  nodes = [na,nb,nc,nd,ne,nf,ng]; edges = [e1,e2,e3,e4,e5,e6,e7,e8]
  graph_string = graph_string + run(nodes,edges,4,8)
  #
  na = Node('A'); nb = Node('B'); nc = Node('C'); nd = Node('D'); ne = Node('E')  ; nf = Node('F'); ng = Node('G')
  e1=Edge(na,nb);e2=Edge(nb,nc);e3=Edge(na,nf,3);e4=Edge(nb,nd);e5=Edge(nd,nf);e6=Edge(nc,ne);
  e7=Edge(nf,ng,2);e8=Edge(ne,ng)
  nodes = [na,nb,nc,nd,ne,nf,ng]; edges = [e1,e2,e3,e4,e5,e6,e7,e8]
  graph_string = graph_string + run(nodes,edges,4,9)
    
  ## Ex 5
  na=Node('A');nb=Node('B');nc=Node('C')
  e1=Edge(na,nb);e2=Edge(nb,na,1);e3=Edge(nc,na)
  nodes=[na,nb,nc];edges=[e1,e2,e3]
  graph_string = graph_string + run(nodes,edges,7,10)
    
  ## Ex 6
  na=Node('A');nb=Node('B');nx=Node('X');ns=Node('S');nd=Node('D');nz=Node('Z')
  e1=Edge(na,nx);e2=Edge(nb,nx);e3=Edge(nx,ns);e4=Edge(nx,nd,1)
  nodes=[na,nb,nd,nx,ns,nz];edges=[e1,e2,e3,e4]
  graph_string = graph_string + run(nodes,edges,4,11)
  
  fout = open('graphs.dot','w')
  fout.write(graph_string)
  fout.close()

