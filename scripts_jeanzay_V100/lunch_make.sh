#!/bin/bash

## JOB INFO
## JOB INFO
#SBATCH --job-name=make_pdi_code
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

## NODE CONFIGURATION
#SBATCH -C v100-16g                 # decommenter pour reserver uniquement des GPU V100 16 Go, quadri GPU
##SBATCH -C v100-32g                 # decommenter pour reserver uniquement des GPU V100 32 Go, quadri GPU
##SBATCH --partition=gpu_p2s         # decommenter pour la partition gpu_p2 (GPU V100 32 Go), octo GPU
##SBATCH --partition=gpu_p2l         # decommenter pour la partition gpu_p2 (GPU V100 32 Go), octo GPU
##SBATCH -C a100                      # 8GPU per node
#SBATCH --qos=qos_gpu-dev # for V100
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10           # nombre de CPU par tache (1/4 du noeud 4-GPU V100 ici)
##SBATCH --cpus-per-task=3           # nombre de CPU par tache pour gpu_p2 (1/8 du noeud 8-GPU V100) partition gpu_p2
##SBATCH --cpus-per-task=8           # nombre de CPU par tache pour gpu_p5 (1/8 du noeud 8-GPU A100)
##SBATCH --cpus-per-task=24          # nombre de CPU par tache pour gpu_p6 (1/4 du noeud 4-GPU H100)
#SBATCH --hint=nomultithread

## JOB ACCOUNTABILITY
#SBATCH -A jyd@v100
#SBATCH --time=00:20:00


./build_all.sh
