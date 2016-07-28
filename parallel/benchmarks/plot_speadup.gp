set terminal png nocrop size 800,600 enhanced font "Verdana, 10" 
set auto x
set auto y
#set logscale y
set grid 
set colorsequence podo
set style data line
set xtic rotate by -0 scale 0
set output 'lanczos.png'
set xlabel 'Number of Vertices'
set ylabel 'Time Spent (s)'
set title 'Figure. Speedup - Parallel Lanczos'
plot for [COL=2:6] 'lanczos.dat' using (column(COL)/$2):xtic(1) ti col lw 2