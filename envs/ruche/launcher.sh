#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=%x_%j.out 
#SBATCH --time=00:10:00 
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=20
#SBATCH --threads-per-core=1
#SBATCH --partition=cpu_short

# All paths are relative to WORKING_DIRECTORY

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
BASE_DIR=${HOME}/bench-in-situ
WORKING_DIR=${BASE_DIR}/working_dir

SCHEFILE=scheduler.json
PREFIX=bench_insitu
DASK_WORKER_NODES=1
SIM_NODES=$(($SLURM_NNODES-2-$DASK_WORKER_NODES))
SIM_PROC=$SIM_NODES

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores

echo "SLURM_NNODES=$SLURM_NNODES"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "SIM_NODES=$SIM_NODES"

# this file must be accessible from every slurm node (i.e.: shared network drive)
source ${BASE_DIR}/envs/ruche/modules.env

# move to working directory 
cd ${WORKING_DIR}

# activate python environment
source deisa/bin/activate

# PDI
source pdi/share/pdi/env.sh 

# dask scheduler
srun -N 1 -n 1 -c 1 -r 0 dask scheduler --scheduler-file=${SCHEFILE} >> ${PREFIX}_dask-scheduler.o &

# Wait for the SCHEFILE to be created
while ! [ -f ${SCHEFILE} ]; do
  sleep 3
done

# dask workers
srun -N ${DASK_WORKER_NODES} -n ${DASK_WORKER_NODES} -c 1 -r 1 dask worker --local-directory /tmp --scheduler-file=${SCHEFILE} >> ${PREFIX}_dask-worker.o &

# insitu
srun -N 1 -n 1 -c 1 -r $(($DASK_WORKER_NODES+1)) python -O in-situ/fft_updated.py >> ${PREFIX}_client.o &
client_pid=$!

# simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} -r $(($DASK_WORKER_NODES+2)) build/main ${BASE_DIR}/envs/ruche/setup.ini ${BASE_DIR}/envs/ruche/io_deisa.yml --kokkos-map-device-id-by=mpi_rank &
simu_pid=$!
wait $simu_pid

