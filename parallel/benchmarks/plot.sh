#!/bin/bash

i=3
for f in lanczos.dat tqli.dat partition.dat
do
	echo "Procs 1024 2048 5120 10240 102400" > $f
	#awk -F" " '{print $2; print "\n";}' ./times/t_16.dat | xargs > $f
	for t in $(ls ./times | sort -n)
	do
	echo $t
	awk -F" " -v i=$i 'FNR == 1 {print $1}{print $i; print "\n";}' ./times/$t | xargs >> $f
	done
	let i++
done
gnuplot < plot.gp
