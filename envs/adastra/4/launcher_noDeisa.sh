#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=res4N_%x_%j.out 
#SBATCH --time=24:00:00 
#SBATCH --nodes=1
#SBATCH --account=cad14985 
#SBATCH --constraint=MI250
#SBATCH --ntasks-per-node=4

# All paths are relative to WORKING_DIRECTORY
SIMU_SIZE=4
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
BASE_DIR=${HOME}/bench-in-situ
WORKING_DIR=${BASE_DIR}/working_dir

PREFIX=bench_insitu
SIM_NODES=1
SIM_PROC=4

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores

echo "SLURM_NNODES=$SLURM_NNODES"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "SIM_NODES=$SIM_NODES"

# this file must be accessible from every slurm node (i.e.: shared network drive)
source ${BASE_DIR}/envs/adastra/modules.env

# set result file path
mkdir -p $SNAPSHOT_FILE_PATH/$SIMU_SIZE
sed -i "s|^prefix=.*|prefix=$SNAPSHOT_FILE_PATH/$SIMU_SIZE/Checkpoint|" ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini

# move to working directory 
cd ${WORKING_DIR}

# PDI
source pdi/share/pdi/env.sh 

# simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} build/main ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini ${BASE_DIR}/envs/adastra/io_chkpt.yml --kokkos-map-device-id-by=mpi_rank &
simu_pid=$!
wait $simu_pid

rm ${SNAPSHOT_FILE_PATH}/${SIMU_SIZE}/*.h5 && rm ${SNAPSHOT_FILE_PATH}/${SIMU_SIZE}/*.xmf