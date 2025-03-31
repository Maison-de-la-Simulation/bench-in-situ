#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=res%%SIMU_SIZEN_%x_%j.out
#SBATCH --time=02:00:00
##################################################
#SBATCH -C partition
##################################################
##SBATCH --ntasks=%%NTASKS, one for each GPU
##SBATCH --nodes=%%SIMU_SIZE/gpus-per-node
#SBATCH --ntasks-per-node=%%NTASKS_PER_NODES
#SBATCH --gres=gpu:%%GPU_PER_NODE
##################################################
#SBATCH --cpus-per-task=%%cpus-per-task
##################################################
#SBATCH --hint=nomultithread
#SBATCH -A %%${IDRPROJ}@${ARCH}

GPU_ARCH=%%GPU_ARCH

# All paths are relative to the WORKING_DIRECTORY of this script
SIMU_SIZE=%%NB_GPU
FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")

##================================================
## Creation of job directory
JOB_GENERATED_DIR=${PWD} ## = ${ENVS_JEANZAY_DIR}/strong_OR_weak_OR_deisa/result_dir/${FORMATED_CUBE_SIZE}/nb_gpu_${FORMATED_SIMU_SIZE}

SIM_NODES=${SLURM_NNODES}
SIM_PROC=${SLURM_NTASKS}

## WORKING DIRECTORY = DIRECTORY WHERE ARE THE DIFFERENT BUILD
# LIBRARY AND EXECUTABLE VARIABLE
PDI_SOURCE=${JOB_GENERATED_DIR}/../../../../working_dir_${GPU_ARCH}/pdi/install/share/pdi/env.sh
MAIN_SIMULATION=${JOB_GENERATED_DIR}/../../../../working_dir_${GPU_ARCH}/sim/build/main

# INPUT FILE VARIABLE
YAML_FILE=%%YAML_FILE

echo "== INFO DIRECTORY"
echo "JOB_GENERATED_DIR=$JOB_GENERATED_DIR"
echo "PDI_SOURCE=$PDI_SOURCE"
echo "MAIN_SIMULATION=$MAIN_SIMULATION"
echo " "

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

echo "== INFO SLURM"
echo "SLURM_NNODES=$SLURM_NNODES"
echo "SLURM_NTASKS=$SLURM_NTASKS"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "SIMU_SIZE=$SIMU_SIZE"
echo " "

# Modules files must be accessible from every slurm node (i.e.: shared network drive)
source ${JOB_GENERATED_DIR}/modules_${GPU_ARCH}.env

# Set result file path
mkdir -p $SNAPSHOT_FILE_PATH/$FORMATED_SIMU_SIZE
sed -i "s|^prefix=.*|prefix=$SNAPSHOT_FILE_PATH/$FORMATED_SIMU_SIZE/Checkpoint|" ${JOB_GENERATED_DIR}/setup.ini

# Move to working directory
cd ${JOB_GENERATED_DIR}

# PDI
source ${PDI_SOURCE}

# Simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} ${MAIN_SIMULATION} ${JOB_GENERATED_DIR}/setup.ini ${JOB_GENERATED_DIR}/${YAML_FILE} --kokkos-map-device-id-by=mpi_rank &
simu_pid=$!
wait $simu_pid

rm ${SNAPSHOT_FILE_PATH}/${FORMATED_SIMU_SIZE}/*.h5 && rm ${SNAPSHOT_FILE_PATH}/${FORMATED_SIMU_SIZE}/*.xmf
