#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=res%%NB_GPU_%x_%j.out
#SBATCH --time=01:00:00
##################################################
#SBATCH -C partition
##################################################
##SBATCH --ntasks=%%NTASKS
##SBATCH --nodes=%%NODES
#SBATCH --ntasks-per-node=%%NTASKS_PER_NODES
#SBATCH --gres=gpu:%%GPU_PER_NODE
##################################################
#SBATCH --cpus-per-task=658
##################################################
#SBATCH --hint=nomultithread
#SBATCH -A %%${IDRPROJ}@${ARCH}

GPU_ARCH=%%GPU_ARCH

# All paths are relative to the DIRECTORY of this script
SIMU_SIZE=%%NB_GPU
FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")

##================================================
## definition of direction
LAUNCH_DIR=${PWD} ## = ${ENVS_JEANZAY_DIR}/strong_OR_weak_OR_deisa/result_dir/${FORMATED_CUBE_SIZE}/nb_gpu_${FORMATED_SIMU_SIZE}

ENVS_JEANZAY_DIR=${LAUNCH_DIR}/../../../.. ## = ${BENCH_ROOT_DIR}
## RESULT DIRECTORY = DIRECTORY WHERE INPUT FILE OF EXECUTABLE

## WORKING DIRECTORY = DIRECTORY WHERE ARE THE DIFFERENT BUILD
# LIBRARY AND EXECUTABLE VARIABLE
PDI_INSTALL_DIR=${ENVS_JEANZAY_DIR}/working_dir_${GPU_ARCH}/pdi/install
MAIN_EXE_DIR=${ENVS_JEANZAY_DIR}/working_dir_${GPU_ARCH}/sim/build

# INPUT FILE VARIABLE
YAML_FILE=%%""

echo "== INFO DIRECTORY"
echo "BASE_DIR=$BASE_DIR"
echo "LAUNCH_DIR=$LAUNCH_DIR"
echo "PDI_INSTALL_DIR=$PDI_INSTALL_DIR"
echo "MAIN_EXE_DIR=$MAIN_EXE_DIR"
echo " "

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

SIM_NODES=${SLURM_NNODES}
SIM_PROC=${SLURM_NTASKS}

echo "== INFO SLURM"
echo "Nb Node(s)=$SLURM_NNODES"
echo "Nb GPU(s) =$SLURM_NTASKS"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "SIMU_SIZE=$SIMU_SIZE"
echo " "

# this file must be accessible from every slurm node (i.e.: shared network drive)
source ${LAUNCH_DIR}/modules_${GPU_ARCH}.env

# set result file path
mkdir -p $SNAPSHOT_FILE_PATH/$SIMU_SIZE
sed -i "s|^prefix=.*|prefix=$SNAPSHOT_FILE_PATH/$SIMU_SIZE/Checkpoint|" ${LAUNCH_DIR}/setup.ini

# move to working directory
cd ${LAUNCH_DIR}

# PDI
source ${PDI_INSTALL_DIR}/share/pdi/env.sh

# simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} ${MAIN_EXE_DIR}/main ${LAUNCH_DIR}/setup.ini ${LAUNCH_DIR}/${YAML_FILE} --kokkos-map-device-id-by=mpi_rank &
simu_pid=$!
wait $simu_pid

rm ${SNAPSHOT_FILE_PATH}/${SIMU_SIZE}/*.h5 && rm ${SNAPSHOT_FILE_PATH}/${SIMU_SIZE}/*.xmf
