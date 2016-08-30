set terminal postscript enhanced color
set auto x
set auto y
set logscale y
set grid 
set style data linespoints
set xtic rotate by -0 scale 0
set xlabel 'Number of processors'
set ylabel 'Time in seconds'

set output 'lanczos.eps'
#set title 'Figure. Elapsed Time - Parallel Lanczos'
plot for [COL=2:7] 'lanczos.dat' using COL:xtic(1) ti col lw 2

set output 'tqli.eps'
#set title 'Figure. Elapsed Time - TQLI'
plot for [COL=2:7] 'tqli.dat' using COL:xtic(1) ti col lw 2

set output 'partition.eps'
#set title 'Figure. Elapsed Time - Parallel Partition'
plot for [COL=2:7] 'partition.dat' using COL:xtic(1) ti col lw 3

