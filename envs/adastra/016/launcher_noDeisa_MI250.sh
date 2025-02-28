#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=res16N_%x_%j.out 
#SBATCH --time=01:00:00 
#SBATCH --nodes=2
#SBATCH --account=cad14985 
#SBATCH --constraint=MI250
##SBATCH --constraint=GENOA
#SBATCH --exclusive
##SBATCH --ntasks-per-node=8
##SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=16
##SBATCH --threads-per-core=1

export MPICH_GPU_SUPPORT_ENABLED=1

# All paths are relative to WORKING_DIRECTORY
SIMU_SIZE=16
BASE_DIR=${PWD}/../..
WORKING_DIR=${BASE_DIR}/working_dir
FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")

PREFIX=bench_insitu
SIM_NODES=2
SIM_PROC=16

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

echo "SLURM_NNODES=$SLURM_NNODES"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "SIM_NODES=$SIM_NODES"

# this file must be accessible from every slurm node (i.e.: shared network drive)
source /lus/home/CT6/cad14985/SHARED/modulesMI250.env

# set result file path
mkdir -p $SNAPSHOT_FILE_PATH/$FORMATED_SIMU_SIZE
sed -i "s|^prefix=.*|prefix=$SNAPSHOT_FILE_PATH/$FORMATED_SIMU_SIZE/Checkpoint|" ${BASE_DIR}/envs/adastra/${FORMATED_SIMU_SIZE}/setup.ini

# move to working directory 
cd ${WORKING_DIR}

# PDI
source ${BASE_DIR}/lib/pdi/build/staging/share/pdi/env.sh

# simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} ${BASE_DIR}/build/main ${BASE_DIR}/envs/adastra/${FORMATED_SIMU_SIZE}/setup.ini ${BASE_DIR}/envs/adastra/io_chkpt.yml --kokkos-map-device-id-by=mpi_rank &
simu_pid=$!
wait $simu_pid

rm ${SNAPSHOT_FILE_PATH}/${FORMATED_SIMU_SIZE}/*.h5 && rm ${SNAPSHOT_FILE_PATH}/${FORMATED_SIMU_SIZE}/*.xmf
