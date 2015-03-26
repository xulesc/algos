#!/usr/bin/python

import numpy as np

CTHRESH = 5
CWINDOWLEN = 15
CWINDOWMID = CWINDOWLEN / 2

def words2dict(fname):
  print 'doing words to dictionary'
  ret = {}
  for line in open(fname):
    for word in line.replace('\n','').split(' '):
      v = ret.get(word)
      if v == None:
	v = 0
      ret[word] = v + 1
  return ret

def writedict(dictionary, fname):
  print 'writing dictionary'
  of = open(fname, 'w'); count = 0
  for v in dictionary:
    of.write('%s %d %d\n' %(v[0], count, v[1]))
    count += 1
  of.close()
  
def readdict(fname):
  print 'reading dictionary'
  ret = {}
  for line in open(fname):
    d = line.replace('\n','').split(' ')
    ret[d[0]] = [int(d[1]), int(d[2])]
  return ret

def makebasemodel(fname, dictionary):
  print 'making base model'
  context = []; dlen = len(dictionary);
  cooc = {}
  for line in open(fname):
    for word in line.replace('\n','').split(' '):
      if word not in dictionary:
	continue
      if len(context) < CWINDOWLEN:
	context.append(dictionary[word][0])
	continue
      main_word_index = context[CWINDOWMID]
      for windex in context:
	if cooc.get(main_word_index) == None:
	  cooc[main_word_index] = {}
	v = cooc[main_word_index]
	if v.get(windex) == None:
	  v[windex] = 0
	cooc[main_word_index][windex] += 1
      context = context[1:]
      context.append(dictionary.get(word)[0])
  return cooc

def writecooc(fname, cooc):
  print 'writing cooc'
  of = open(fname, 'w')
  for k, v in cooc.iteritems():
    for p, c in v.iteritems():
      of.write('%d %d %d\n' %(k, p, c))
  of.close()
  
def readcooc(fname):
  print 'reading cooc'
  cooc = {}
  for line in open(fname):
    (w1,w2,c) = line.replace('\n','').split(' ')
    w1 = int(w1); w2 = int(w2); c = int(c)
    if cooc.get(w1) == None:
      cooc[w1] = {}
    cooc[w1][w2] = c
  return cooc	  
  
def unit_vector(vector):
  return vector / np.linalg.norm(vector)

def cos_distance(v1, v2):
  v1_u = unit_vector(v1)
  v2_u = unit_vector(v2)
  return np.dot(v1_u, v2_u)  

def angle_between(v1, v2):
  angle = np.arccos(cos_distance(v1, v2))
  if np.isnan(angle):
    if (v1_u == v2_u).all():
      return 0.0
    else:
      return np.pi
  return angle

def compare(w1_np_arr, w2, dictionary, cooc):
  w2_index = dictionary[w2][0]
  w2_vec = cooc[w2_index]
  dlen = len(dictionary)
  rmssq = 0
  w2_np_arr = np.zeros(dlen)
  for k, v in w2_vec.iteritems():
    w2_np_arr[k] = v
  return cos_distance(w1_np_arr, w2_np_arr)
    
###################### PROCESSING ##############################################
infile_name = 'text8'
dictionary_name = '%s.dictionary' %infile_name
base_model_name = '%s.base.model' %infile_name
## words to dict
dictionary = words2dict(infile_name)
dictionary = filter(lambda x : x[1] >= CTHRESH, dictionary.iteritems() )
writedict(dictionary, dictionary_name)
## make base model
dictionary = readdict(dictionary_name)
cooc = makebasemodel(infile_name, dictionary)
writecooc(base_model_name, cooc)
## find nearest
dictionary = readdict(dictionary_name)
cooc = readcooc(base_model_name)
of = open('run.log', 'w')
test_word_str = 'schizophrenia'
##
w1_index = dictionary[test_word_str][0]
w1_vec = cooc[w1_index]
dlen = len(dictionary)
w1_np_arr = np.zeros(dlen)
for k, v in w1_vec.iteritems():
  w1_np_arr[k] = v
##
print 'comparing'
for k, v in dictionary.iteritems():
  dist = compare(w1_np_arr,k,dictionary,cooc)
  outstr = '>%s %s %f' %(test_word_str, k, dist)
  of.write('%s\n' %outstr)
of.close()
##
#print 'france spain %f' %compare('france','spain',dictionary,cooc)
#print 'france italy %f' %compare('france','italy',dictionary,cooc)
#print 'france french %f' %compare('france','french',dictionary,cooc)
#print 'france german %f' %compare('france','german',dictionary,cooc)
#print 'france belgium %f' %compare('france','belgium',dictionary,cooc)
#print 'france commune %f' %compare('france','commune',dictionary,cooc)
#print 'france provence %f' %compare('france','provence',dictionary,cooc)
#print 'france netherlands %f' %compare('france','netherlands',dictionary,cooc)
#print 'france alsace %f' %compare('france','alsace',dictionary,cooc)
#print 'france paris %f' %compare('france','paris',dictionary,cooc)
##
#print 'india pakistan %f' %compare('india','pakistan',dictionary,cooc)
#print 'india hyderabad %f' %compare('india','hyderabad',dictionary,cooc)
#print 'india assam %f' %compare('india','assam',dictionary,cooc)
#print 'india cochin %f' %compare('india','cochin',dictionary,cooc)
#print 'india pune %f' %compare('india','pune',dictionary,cooc)
#
#print 'food foods %f' %compare('food','foods',dictionary,cooc)
#print 'food meat %f' %compare('food','meat',dictionary,cooc)
#print 'food nutrition %f' %compare('food','nutrition',dictionary,cooc)
#print 'food undigested %f' %compare('food','undigested',dictionary,cooc)
#print 'food shellfish %f' %compare('food','shellfish',dictionary,cooc)
