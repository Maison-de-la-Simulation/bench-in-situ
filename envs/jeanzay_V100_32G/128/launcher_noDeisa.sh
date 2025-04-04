#!/bin/bash

#SBATCH --job-name=bench_insitu
#SBATCH --output=res128%x_%j.out
#SBATCH --time=01:00:00
##################################################
##SBATCH -C v100-16g                  # decommenter pour reserver uniquement des GPU V100 16 Go, quadri GPU
#SBATCH -C v100-32g                 # decommenter pour reserver uniquement des GPU V100 32 Go, quadri GPU
##SBATCH --partition=gpu_p2s         # decommenter pour la partition gpu_p2 (GPU V100 32 Go), octo GPU
##SBATCH --partition=gpu_p2l         # decommenter pour la partition gpu_p2 (GPU V100 32 Go), octo GPU
##SBATCH -C a100                     # 8GPU per node
##################################################
#SBATCH --ntasks=128
#SBATCH --nodes=32                    ##perhaps need to comment ?
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
##################################################
#SBATCH --cpus-per-task=10           # nombre de CPU par tache (1/4 du noeud 4-GPU V100 ici)
##SBATCH --cpus-per-task=3           # nombre de CPU par tache pour gpu_p2 (1/8 du noeud 8-GPU V100) partition gpu_p2
##SBATCH --cpus-per-task=8           # nombre de CPU par tache pour gpu_p5 (1/8 du noeud 8-GPU A100)
##SBATCH --cpus-per-task=24          # nombre de CPU par tache pour gpu_p6 (1/4 du noeud 4-GPU H100)
##################################################
#SBATCH --hint=nomultithread
#SBATCH -A jyd@v100

# All paths are relative to WORKING_DIRECTORY

NODES_ARCH_COMPILATION="_V100"
PDI_MHD_NODES_ARCH="_V100_32G"

# All paths are relative to WORKING_DIRECTORY
SIMU_SIZE=${SLURM_NTASKS}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
BASE_DIR=${WORK}/numpex/bench-in-situ
WORKING_DIR=${BASE_DIR}/working_dir${NODES_ARCH_COMPILATION}

PREFIX=bench_insitu
SIM_NODES=${SLURM_NNODES}
SIM_PROC=${SLURM_NTASKS}

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
#export OMP_PROC_BIND=spread
#export OMP_PLACES=threads


echo "SLURM_NNODES=$SLURM_NNODES"
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "SIM_NODES=$SIM_NODES"

# this file must be accessible from every slurm node (i.e.: shared network drive)
source ${BASE_DIR}/envs/jeanzay${PDI_MHD_NODES_ARCH}/modules.env

# set result file path
mkdir -p $SNAPSHOT_FILE_PATH/$SIMU_SIZE
sed -i "s|^prefix=.*|prefix=$SNAPSHOT_FILE_PATH/$SIMU_SIZE/Checkpoint|" ${BASE_DIR}/envs/jeanzay${PDI_MHD_NODES_ARCH}/${SIMU_SIZE}/setup.ini

# move to working directory
cd ${WORKING_DIR}

# PDI
source pdi/share/pdi/env.sh

# simulation
srun -N ${SIM_NODES} -n ${SIM_PROC} build/main ${BASE_DIR}/envs/jeanzay${PDI_MHD_NODES_ARCH}/${SIMU_SIZE}/setup.ini ${BASE_DIR}/envs/jeanzay${PDI_MHD_NODES_ARCH}/io_chkpt.yml --kokkos-map-device-id-by=mpi_rank &
simu_pid=$!
wait $simu_pid

rm ${SNAPSHOT_FILE_PATH}/${SIMU_SIZE}/*.h5 && rm ${SNAPSHOT_FILE_PATH}/${SIMU_SIZE}/*.xmf
