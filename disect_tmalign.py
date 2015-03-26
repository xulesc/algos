#!/usr/bin/python

import numpy as np

def _get_dom_name(string):
  return string[string.rindex('/') + 1 :]
  
def _read_alignment(string1, string2):
  loc1=0;loc2=0;edges=[]
  for e1,e2 in zip(string1, string2):
    if e1 != '-' and e2 != '-':
      edges.append([loc1,loc2])
      loc1 += 1; loc2 += 1
    elif e1 == '-':
      loc2 += 1
    elif e2 == '-':
      loc1 += 1
  return edges

def _get_rmat_row(string):
  data =  string.split(' ')
  found = 0
  for d in data:
    if len(d) == 0:
      continue
    found += 1
    if found == 1:
      continue
    if found == 2:
      t = float(d)
    elif found == 3:
      u0 = float(d)
    elif found == 4:
      u1 = float(d)
    elif found == 5:
      u2 = float(d)
  return [t, u0, u1, u2]

def read_tmalign_output(fname):
  rmat = []; seq_ready = False; cnt = 0; l1 = ''; l2 = ''
  for line in open(fname, 'r'):
    line = line.replace('\n','')
    if line.startswith('Name of Chain_1'):
      dom1 = _get_dom_name(line)
    if line.startswith('Name of Chain_2'):
      dom2 = _get_dom_name(line)
    if line.startswith(' 1'):
      rmat.append(_get_rmat_row(line))
    if line.startswith(' 2'):
      rmat.append(_get_rmat_row(line))
    if line.startswith(' 3'):
      rmat.append(_get_rmat_row(line))
    if len(line) > 3 and line[2] == ':': #.startswith("\(\":\""):
      seq_ready = True
      continue
    cnt += (seq_ready == True)
    if seq_ready == True and cnt == 1:
      l1 = line
    if seq_ready == True and cnt == 3:
      l2 = line
  print l1
  print l2    
  return [dom1, dom2, np.array(rmat), _read_alignment(l1,l2)]

if __name__ == "__main__":
  print read_tmalign_output('delme.1aa9.1ash.out')
  print 'done'