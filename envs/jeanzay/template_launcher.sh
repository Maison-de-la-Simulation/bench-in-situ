#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=res%%SIMU_SIZE_N_%x_%j.out
#SBATCH --time=02:00:00
##################################################
#SBATCH -C partition
##################################################
##SBATCH --ntasks=%%NTASKS #NB_GPU
##SBATCH --nodes=%%NODES   #NB_NODES
#SBATCH --ntasks-per-node=%%NTASKS_PER_NODES
#SBATCH --gres=gpu:%%GPU_PER_NODE
##################################################
#SBATCH --cpus-per-task=%%CPUS-PER-TASKS
##################################################
#SBATCH --hint=nomultithread
#SBATCH -A %%${IDRPROJ}@${ARCH}

GPU_ARCH=%%GPU_ARCH

##================================================
## Directory and filename
##================================================
# All paths are relative to the working directory "JOB_GENERATED_DIR" of this script
JOB_GENERATED_DIR=${PWD} ## = ${ENVS_JEANZAY_DIR}/strong_OR_weak_OR_deisa/result_dirname/${FORMATED_CUBE_SIZE}/nb_gpu_${FORMATED_SIMU_SIZE}

# LIBRARY AND EXECUTABLE VARIABLE
PDI_ENV_SCRIPT=%%PDI_INSTALL_DIR%%/share/pdi/env.sh
MAIN_SIMULATION=%%MAIN_EXE_DIR%%/main

# INPUT FILE VARIABLE
YAML_FILE=%%YAML_FILE

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

SIM_NODES=${SLURM_NNODES}
SIM_PROC=${SLURM_NTASKS}

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
source ${PDI_ENV_SCRIPT}

# Simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} ${MAIN_SIMULATION} ${JOB_GENERATED_DIR}/setup.ini ${JOB_GENERATED_DIR}/${YAML_FILE} --kokkos-map-device-id-by=mpi_rank &
simu_pid=$!
wait $simu_pid

rm ${SNAPSHOT_FILE_PATH}/${FORMATED_SIMU_SIZE}/*.h5 && rm ${SNAPSHOT_FILE_PATH}/${FORMATED_SIMU_SIZE}/*.xmf
rm -rf ${SNAPSHOT_FILE_PATH}/${FORMATED_SIMU_SIZE}