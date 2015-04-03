#!/usr/bin/python

import sys

inf=open(sys.argv[1])
ouf=open(sys.argv[2],'w')

##
similarities = {}
for line in inf:
  data = line.replace('\n','').split(' ')
  if similarities.get(data[0]) == None:
    similarities[data[0]] = {}
  similarities[data[0]][data[1]] = float(data[2])
##
for k, v in similarities.iteritems():
  minv = min(v.values())
  maxv = max(v.values())
  orng = maxv - minv
  for k1, v1 in v.iteritems():
    ouf.write('%s %s %f\n' %(k, k1, ((v1 - minv) / orng)))
##

inf.close()
ouf.close()

##