#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=%x_%j.out
#SBATCH --time=00:10:00
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=4
#SBATCH --nodes=4
#SBATCH -C a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH -A wuc@a100

# All paths are relative to WORKING_DIRECTORY

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
BASE_DIR=${WORK}/numpex/bench-in-situ
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
source ${BASE_DIR}/envs/jeanzay/modules.env

# move to working directory
cd ${WORKING_DIR}

# activate python environment
source deisa/bin/activate

# PDI
source pdi/share/pdi/env.sh

export LD_LIBRARY_PATH="/gpfslocalsup/pub/anaconda-py3/2023.09/envs/python-3.11.5/lib:$LD_LIBRARY_PATH"
# dask scheduler
srun -N 1 -n 1 -c 1 -r 0 dask scheduler --protocol tcp --scheduler-file=${SCHEFILE} >> ${PREFIX}_dask-scheduler.o &

# Wait for the SCHEFILE to be created
while ! [ -f ${SCHEFILE} ]; do
  sleep 3
done

# dask workers
srun -N ${DASK_WORKER_NODES} -n ${DASK_WORKER_NODES} -c 1 -r 1 dask worker --protocol tcp --local-directory /tmp --scheduler-file=${SCHEFILE} >> ${PREFIX}_dask-worker.o &

# insitu
srun -N 1 -n 1 -c 1 -r $(($DASK_WORKER_NODES+1)) python -O in-situ/fft_updated.py >> ${PREFIX}_client.o &
client_pid=$!

# simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} -r $(($DASK_WORKER_NODES+2)) build/main ${BASE_DIR}/envs/jeanzay/setup.ini ${BASE_DIR}/envs/jeanzay/io_deisa.yml --kokkos-map-device-id-by=mpi_rank &
simu_pid=$!
wait $simu_pid