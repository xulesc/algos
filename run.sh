#!/bin/sh

fs=`ls pdb_data`
for f1 in $fs; do
    for f2 in $fs; do
        d1=`echo $f1 | sed -e 's/.pdb//g'`
        d2=`echo $f2 | sed -e 's/.pdb//g'`
        python test_psalign.py $d1 $d2
    done
done
