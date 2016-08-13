set terminal png nocrop size 800,600 enhanced font "Verdana, 11" 
set auto x
set auto y
#set logscale y
set grid 
set key left top
set colorsequence podo
set style data linespoints
set xtic rotate by -0 scale 0
set xlabel 'Number of processors'
set ylabel 'Speed-up'

set output 'speedup.png'
set title 'Figure. Speed-up - Parallel Partition'
plot for [COL=2:6] 'speedup.dat' using COL:xtic(1) ti col lw 2

