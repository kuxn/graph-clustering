#!/bin/bash

i=3
for f in ./benchmarks/lanczos.dat ./benchmarks/tqli.dat ./benchmarks/partition.dat
do
#paste ./*.dat | awk -F" " '{printf("%s ", $2);for(x=index;x<=NF;x+=5)printf("%s ",$x); printf("\n");}' > $f
paste ./benchmarks/times/*.dat | awk -F" " -v i=$i 'BEGIN{ORS="";}{print $2," "; for(x=i;x<=NF;x+=5)print $x," "; print "\n";}' > $f
let i++
#sed -i '' '1d' $f
sed -i '' '1i\
Vertices  "4 procs" "8 procs"  "16 procs"   "32 procs"  "64 procs"
' $f
done
gnuplot < ./benchmarks/plot.gp
