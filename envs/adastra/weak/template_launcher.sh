#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=res%%SIMU_SIZEN_%x_%j.out
#SBATCH --time=02:00:00
##################################################
#SBATCH --constraint=%%constraint
##SBATCH --constraint=GENOA
##################################################
#SBATCH --nodes=%%SIMU_SIZE/gpus-per-node
#SBATCH --exclusive
##SBATCH --ntasks-per-node=8, one for each GCD
##SBATCH --gpus-per-node=8 GCD for MI250, 4 for MI300
#SBATCH --cpus-per-task=%%cpus-per-task
##SBATCH --threads-per-core=1
##################################################
#SBATCH --account=%%ACTIVE_PROJECT

export MPICH_GPU_SUPPORT_ENABLED=1

##================================================
## Directory and filename
##================================================
# All paths are relative to the working directory "JOB_GENERATED_DIR" of this script
ROOT_DIR=${PWD}/../../..
JOB_GENERATED_DIR=%%JOB_GENERATED_DIR

# LIBRARY AND EXECUTABLE VARIABLE
PDI_ENV_SCRIPT=${ROOT_DIR}/lib/install_pdi/share/pdi/env.sh
MAIN_SIMULATION=${ROOT_DIR}/simulation/build/main

echo "== INFO DIRECTORY AND FILENAME"
echo "JOB_GENERATED_DIR=$JOB_GENERATED_DIR"
echo "PDI_ENV_SCRIPT=$PDI_ENV_SCRIPT"
echo "MAIN_SIMULATION=$MAIN_SIMULATION"
echo " "

##=================================================
## Job parameters
##=================================================
SIMU_SIZE=%%NB_GPU
FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")

SIM_NODES=%%SIM_NODES
SIM_PROC=%%SIMU_SIZE

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

echo "== INFO SLURM"
echo "SLURM_NNODES=$SLURM_NNODES"
# echo "SLURM_NTASKS=$SLURM_NTASKS"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "SIMU_SIZE=$SIMU_SIZE"
echo " "

# Modules files must be accessible from every slurm node (i.e.: shared network drive)
source ${JOB_GENERATED_DIR}/modules.env
source ${JOB_GENERATED_DIR}/modulesGPU%%constraint.env

# Set result file path
mkdir -p $SNAPSHOT_FILE_PATH/$FORMATED_SIMU_SIZE
sed -i "s|^prefix=.*|prefix=$SNAPSHOT_FILE_PATH/$FORMATED_SIMU_SIZE/Checkpoint|" ${JOB_GENERATED_DIR}/setup.ini

# Move to the working directory
cd ${JOB_GENERATED_DIR}

# PDI
source ${PDI_ENV_SCRIPT}

# rocm-smi
##rocm-smi > rocm-monitor${FORMATED_SIMU_SIZE}.csv

# Simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} ${ROOT_DIR}/simulation/build/main ${JOB_GENERATED_DIR}/setup.ini ${JOB_GENERATED_DIR}/io_chkpt.yml --kokkos-map-device-id-by=mpi_rank &
simu_pid=$! &
rocm-smi > rocm-monitor${FORMATED_SIMU_SIZE}.csv
wait $simu_pid

rm ${SNAPSHOT_FILE_PATH}/${FORMATED_SIMU_SIZE}/*.h5 && rm ${SNAPSHOT_FILE_PATH}/${FORMATED_SIMU_SIZE}/*.xmf
rm -rf ${SNAPSHOT_FILE_PATH}/${FORMATED_SIMU_SIZE}
