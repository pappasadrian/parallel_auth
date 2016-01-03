#!/bin/bash
#PBS -­q auth
#PBS ­-N mpi-knn
#PBS ­-j oe
#PBS ­-l nodes=1:ppn=4

cd $PBS_O_WORKDIR

export NP=$(cat $PBS_NODEFILE | wc -l)

export I2G_MPI_TYPE=mpich2
export I2G_MPI_APPLICATION=knn
export I2G_MPI_APPLICATION_ARGS="q c p b"
echo $I2G_MPI_APPLICATION_ARGS
$I2G_MPI_START
