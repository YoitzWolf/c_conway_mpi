#!/bin/bash
#PBS -q smm
#PBS -l nodes=4:ppn=2

LIB_OLD=$LD_LIBRARY_PATH
PATH_OLD=$PATH

cd ~/z-openmpi-new

echo "CUSTOM OPEN MPI"
	export LD_LIBRARY_PATH=$HOME/OPENMPI/lib:$LIB_OLD
	export 	PATH=$HOME/OPENMPI/bin:$PATH_OLD
	make conway
	mpirun -N 2 (printf '100 500 1000 100' | ./conway) # << 2 processes per each node
echo "CUSTOM OPEN MPI DONE"

echo " "
