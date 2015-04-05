#!/usr/bin/python

import sys
import struct

CLEAR_CODE = '|'
END_CODE = '<'
ALPHABET = "0123456789 \n%s%s" %(CLEAR_CODE, END_CODE)
INDEX_BIT_CNT = 8
MAX_DICT_SIZE = 2 ** INDEX_BIT_CNT

def bitstobytes(bits):
    """
    Interprets an indexable list of booleans as bits, MSB first, to be
    packed into a list of integers from 0 to 256, MSB first, with LSBs
    zero-padded. Note this padding behavior means that round-trips of
    bytestobits(bitstobytes(x, width=W)) may not yield what you expect
    them to if W % 8 != 0

    Does *NOT* pack the returned values into a bytearray or the like.

    >>> import lzw
    >>> bitstobytes([0, 0, 0, 0, 0, 0, 0, 0, "Yes, I'm True"]) == [ 0x00, 0x80 ]
    True
    >>> bitstobytes([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0]) == [ 0x01, 0x30 ]
    True
    """
    ret = []
    nextbyte = 0
    nextbit = 7
    for bit in bits:
        if bit:
            nextbyte = nextbyte | (1 << nextbit)

        if nextbit:
            nextbit = nextbit - 1
        else:
            ret.append(nextbyte)
            nextbit = 7
            nextbyte = 0

    if nextbit < 7: ret.append(nextbyte)
    return ret

## clear dictionary
def _reset_dict(): return dict(map(lambda x : (ALPHABET[x], x), xrange(len(ALPHABET))))
## int to bits
def _encode(index): return map(lambda x : int(x), (INDEX_BIT_CNT*"0"+bin(index)[2:])[-INDEX_BIT_CNT:])
## bits to int 
def _decode(bits): return int(bits,2)
  
def encode(bytes):
  buffer = ''
  dictionary = _reset_dict()
  encoded = []; 
  for b in bytes:
    if buffer + b in dictionary:
      buffer = buffer + b
    else:
      encoded.append(_encode(dictionary[buffer]))
      dictionary[buffer + b] = len(dictionary)
      buffer = b
    if len(dictionary) == MAX_DICT_SIZE:
      encoded.append(dictionary[CLEAR_CODE])
      dictionary = _reset_dict()
      buffer = b
  encoded.append(_encode(dictionary[END_CODE]))
  return encoded
  
def _unencoded(encoded): return map(lambda x : _decode(x), encoded)
  
if __name__ == '__main__':
  #print _reset_dict()
  #print _encode(1)
  #print _decode(_encode(2))
  #tststr = "0 2\n1 306\n1 309\n3 7\n3 45\n4 6\n4 7\n4 45\n5 7\n5 8\n5 45\n5 46\n5 47\n6 8\n7 9\n8 10\n8 47\n8 281\n8 282\n9 46\n9 47\n9 281\n9 282\n9 283\n10 46\n10 47"
  tststr = ''.join([x for x in open(sys.argv[1]) if '#' not in x])
  print len(tststr)
  enc = encode(tststr)
  enc = map(lambda x: map(lambda y : struct.pack("B",y), bitstobytes(x)), enc)
  print(len(enc))


