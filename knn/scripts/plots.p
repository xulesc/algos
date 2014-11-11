set term postscript eps enhanced color
set yrange [0:]

set datafile separator ","

set key font ",11"
set xtics font ",10"
set ytics font ",10"
set key left top

set out "plot_data/cc.eps"
set ylabel "Speedup" font ",13"
set xlabel "Number of threads" font ",13"
stats 'plot_data/non-kdtree-cc.dat' u 1:2
max_y = STATS_max_y
set xrange [0:STATS_max_x + 1]
set arrow 1 from 4,0 to 4,12 nohead lt 3
set label 1 "Number of actual processors" at 3.8,0.7 rotate by 90 font ",11"  
set arrow 2 from 8,0 to 8,12 nohead lt 4
set label 2 "Number of actual cores" at 7.8,4 rotate by 90 font ",11"  
plot 'plot_data/kdtree-cc.dat' u 1:(max_y/$2) with line title 'kdtree-KNN', \
     'plot_data/non-kdtree-cc.dat' u 1:(max_y/$2) with line title 'simple-KNN';


set out "plot_data/k.eps"
set ylabel "Slow down factor" font ",13"
set xlabel "Number of neighbors considered (k)" font ",13"
unset arrow 1
unset arrow 2
unset label 1
unset label 2
stats 'plot_data/kdtree-k.dat' u 1:2
max_y = STATS_max_y
set xrange [0:STATS_max_x + 1]
plot 'plot_data/kdtree-k.dat' u 1:(max_y,$2) with line title '';

##
