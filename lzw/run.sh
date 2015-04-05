#!/bin/sh

DATA_DIR="/media/lvm/xule/workspace/noc_mcpsc/data/contact_maps_chew_kedem/"
SIM_FILE="similarities"
SIM_FILE_F="$SIM_FILE.scaled"
CMD1="python lzw.py"
CMD2="./lzw"
CMD3="python mylzw.py"
CMD=$CMD2

process()
{
cat $DATA_DIR/$1 > t1
cat $DATA_DIR/$2 > t2
cat t1 > t1t2; cat t2 >> t1t2
cat t2 > t2t1; cat t1 >> t2t1

#x=$((`cat t1 | $CMD | wc -c`))
#y=$((`cat t2 | $CMD | wc -c`))
#xy=$((`cat t1t2 | $CMD | wc -c`))
#yx=$((`cat t2t1 | $CMD | wc -c`))
x=$((`$CMD t1`))
y=$((`$CMD t2`))
xy=$((`$CMD t1t2`))
yx=$((`$CMD t2t1`))

echo "" | awk -v nm="$1 $2" -v x=$x -v y=$y -v xy=$xy -v yx=$yx '{d=(x<y)?y:x;n=(xy-x<yx-y)?yx-y:xy-x}END{print nm,n/d}'
}

pairwise_similarities()
{
rm $SIM_FILE
for f1 in `ls $DATA_DIR/`; do
    for f2 in `ls $DATA_DIR/`; do
        process $f1 $f2 >> $SIM_FILE
    done;
done
}

####
### create pariwise sim
pairwise_similarities
### scale pairwise sim
python scale.py $SIM_FILE $SIM_FILE_F

