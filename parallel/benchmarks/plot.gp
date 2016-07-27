set terminal png nocrop size 800,600 enhanced font "Verdana, 10" 
set auto x
set auto y
#set logscale y
set grid 
set colorsequence podo
set style data line
set xtic rotate by -0 scale 0
set output './benchmarks/lanczos.png'
set xlabel 'Number of Vertices'
set ylabel 'Time Spent (s)'
set title 'Figure. Time Spent - Parallel Lanczos'
plot for [COL=2:6] './benchmarks/lanczos.dat' using COL:xtic(1) ti col lw 2

set output './benchmarks/tqli.png'
set xlabel 'Number of Vertices'
set ylabel 'Time Spent (s)'
set title 'Figure. Time Spent - Parallel TQLI'
plot for [COL=2:6] './benchmarks/tqli.dat' using COL:xtic(1) ti col lw 2

set output './benchmarks/partition.png'
set xlabel 'Number of Vertices'
set ylabel 'Time Spent (s)'
set title 'Figure. Time Spent - Parallel Partition'
plot for [COL=2:6] './benchmarks/partition.dat' using COL:xtic(1) ti col lw 2


