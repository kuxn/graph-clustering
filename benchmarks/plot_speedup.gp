set terminal postscript enhanced color
set auto x
set auto y
set grid 
set key left top
set style data linespoints
set xtic rotate by -0 scale 0
set xlabel 'Number of processors'
set ylabel 'Speed-up'

set output 'speedup.eps'
#set title 'Figure. Speed-up - Parallel Partition'
plot for [COL=2:7] 'speedup.dat' using COL:xtic(1) ti col lw 3

