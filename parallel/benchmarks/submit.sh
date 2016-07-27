#!/bin/bash

mkdir data
mkdir times
for i in 2 4 8 16 32 64
do
sbatchname="sbatch_$i.sh"
filename="$i.dat"
##SBATCH -t $(echo "$i / 4" | bc):00:00
cat > $sbatchname << EOL
#!/bin/bash
#SBATCH -n $i
#SBATCH -p compute
#SBATCH -t 2:00:00
#SBATCH -U mschpc
#SBATCH -J run
#SBATCH --reservation=application

# load the modules required
module unload deprecated old/gcc/64/4.6.3
module load apps libs cports
module load boost/1_58_0_openmpi_1.8.6-gnu
module load gcc/4.8.2-gnu

for v in 1024 2048 5120 10240 102400 1024000
do
    time mpirun -np $i ../main -f ./test/par_test_\$v.dot -v \$v -g -o >> ./data/$filename
done

EOL
	sbatch $sbatchname
done
