#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=A100_res64%x_%j.out
##SBATCH --error=A100_res64%x_%j.err
#SBATCH --time=00:20:00
##################################################
##SBATCH -C v100-16g                 # decommenter pour reserver uniquement des GPU V100 16 Go, quadri GPU
##SBATCH -C v100-32g                 # decommenter pour reserver uniquement des GPU V100 32 Go, quadri GPU
##SBATCH --partition=gpu_p2s         # decommenter pour la partition gpu_p2 (GPU V100 32 Go), octo GPU
##SBATCH --partition=gpu_p2l         # decommenter pour la partition gpu_p2 (GPU V100 32 Go), octo GPU
#SBATCH -C a100                      # 8GPU per node
##################################################
#SBATCH --ntasks=64
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:8
##################################################
##SBATCH --cpus-per-task=10         # nombre de CPU par tache (1/4 du noeud 4-GPU V100 ici)
##SBATCH --cpus-per-task=3          # nombre de CPU par tache pour gpu_p2 (1/8 du noeud 8-GPU V100) partition gpu_p2
#SBATCH --cpus-per-task=8           # nombre de CPU par tache pour gpu_p5 (1/8 du noeud 8-GPU A100)
##SBATCH --cpus-per-task=24         # nombre de CPU par tache pour gpu_p6 (1/4 du noeud 4-GPU H100)
##################################################
#SBATCH --hint=nomultithread
##SBATCH --exclusive
#SBATCH -A jyd@a100

# All paths are relative to WORKING_DIRECTORY

GPU_ARCH="A100"
NODES_ARCH_COMPILATION="_A100"
PDI_MHD_NODES_ARCH=""

# All paths are relative to WORKING_DIRECTORY
SIMU_SIZE=${SLURM_NTASKS}
FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")

##================================================
## definition of direction
BENCH_ROOT_DIR=${PWD}
BASE_DIR=${BENCH_ROOT_DIR}

## WORKING DIRECTORY
WORKING_DIR=${BASE_DIR}/working_dir${NODES_ARCH_COMPILATION}

## DIRECTORY USE TO SOURCE FILE
MODULES_DIR=${BASE_DIR}/envs/${PDI_MHD_NODES_ARCH}
PDI_INSTALL_DIR=${WORKING_DIR}

## DIRECTORY OF INPUT FILE
SETUP_INI_DIR=${BASE_DIR}/envs/${PDI_MHD_NODES_ARCH}/${FORMATED_SIMU_SIZE}
YAML_FILE_DIR=${BASE_DIR}/envs/${PDI_MHD_NODES_ARCH}

## DIRECTORY OF EXECUTABLE OF PHYSICS SIMULATION CODE
MAIN_EXE_DIR=${WORKING_DIR}/build

echo "== INFO DIRECTORY"
echo "BASE_DIR=$BASE_DIR"
echo "MODULES_DIR=$MODULES_DIR"
echo "PDI_INSTALL_DIR=$PDI_INSTALL_DIR"
echo "YAML_FILE_DIR=$YAML_FILE_DIR"
echo "SETUP_INI_DIR=$SETUP_INI_DIR"
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
source ${MODULES_DIR}/modules_${GPU_ARCH}.env

# set result file path
mkdir -p $SNAPSHOT_FILE_PATH/$SIMU_SIZE
sed -i "s|^prefix=.*|prefix=$SNAPSHOT_FILE_PATH/$SIMU_SIZE/Checkpoint|" ${SETUP_INI_DIR}/setup.ini

# move to working directory
cd ${WORKING_DIR}

# PDI
source ${PDI_INSTALL_DIR}/pdi/share/pdi/env.sh

# simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} ${MAIN_EXE_DIR}/main ${SETUP_INI_DIR}/setup.ini ${YAML_FILE_DIR}/io_chkpt.yml --kokkos-map-device-id-by=mpi_rank &
simu_pid=$!
wait $simu_pid

rm ${SNAPSHOT_FILE_PATH}/${SIMU_SIZE}/*.h5 && rm ${SNAPSHOT_FILE_PATH}/${SIMU_SIZE}/*.xmf
