#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=res1_%x_%j.out
#SBATCH --time=01:00:00
##################################################
#SBATCH -C partition
##################################################
##SBATCH --ntasks=1
##SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
##################################################
#SBATCH --cpus-per-task=658
##################################################
#SBATCH --hint=nomultithread
#SBATCH -A ${IDRPROJ}@${IDR_ARCH}

# All paths are relative to WORKING_DIRECTORYIDRPROJ
GPU_ARCH="V100"
NODES_ARCH_COMPILATION="_V100"

# All paths are relative to WORKING_DIRECTORY
SIMU_SIZE=1 #${SLURM_NTASKS}
FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")

##================================================
## definition of direction
LAUNCH_DIR=${PWD} ## = ${BENCH_ROOT_DIR}/result_dir/${FORMATED_SIMU_SIZE}

BENCH_ROOT_DIR=${LAUNCH_DIR}/../.. ## = ${BENCH_ROOT_DIR}
## RESULT DIRECTORY = DIRECTORY WHERE INPUT FILE OF EXECUTABLE

## WORKING DIRECTORY = DIRECTORY WHERE ARE THE DIFFERENT BUILD
# LIBRARY AND EXECUTABLE VARIABLE
PDI_INSTALL_DIR=${BENCH_ROOT_DIR}/working_dir${NODES_ARCH_COMPILATION}
MAIN_EXE_DIR=${BENCH_ROOT_DIR}/working_dir${NODES_ARCH_COMPILATION}/build

# INPUT FILE VARIABLE
YAML_FILE="" ##

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
source ${PDI_INSTALL_DIR}/pdi/share/pdi/env.sh

# simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} ${MAIN_EXE_DIR}/main ${LAUNCH_DIR}/setup.ini ${LAUNCH_DIR}/${YAML_FILE} --kokkos-map-device-id-by=mpi_rank &
simu_pid=$!
wait $simu_pid

rm ${SNAPSHOT_FILE_PATH}/${SIMU_SIZE}/*.h5 && rm ${SNAPSHOT_FILE_PATH}/${SIMU_SIZE}/*.xmf
