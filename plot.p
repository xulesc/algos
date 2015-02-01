#
set term png
set output "sq.bench.png"
set key right bottom
set title "SVD vs QCP for superimposing structures" font ",14"
set xlabel "Length of residue" font ",14"
set ylabel "Time (msec)" font ",14"
plot 'sq.bench.out.log' using 1:2 title 'SVD' linecolor rgb "black", 'sq.bench.out.log' using 1:3 title 'QCP' linecolor rgb "black"
