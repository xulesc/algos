#!/bin/sh

grep '#kdtree-cc' log_bench | awk '{print $5 $17}' > plot_data/kdtree-cc.dat
grep '#non-kdtree-cc' log_bench | awk '{print $5 $17}' > plot_data/non-kdtree-cc.dat
grep '#kdtree-k' log_bench | awk '{print $9 $17}' > plot_data/kdtree-k.dat

##

