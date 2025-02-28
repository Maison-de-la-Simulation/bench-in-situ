#!/bin/bash
#SBATCH --account=cad14985
#sbatch --output=1NSizeBenh.out
#SBATCH --job-name=deisaBenchmark
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=12:00:00

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
BASE_DIR=${HOME}/bench-in-situ
WORKING_DIR=${BASE_DIR}/working_dir
SIMU_SIZE=1

cd ${WORKING_DIR}
rm *.h5
rm *.xmf
cd ..


for  ((cubeSize=64; i<=2048; i+=128)) ; do
    sed -i "s/^nx=[0-9]*$/nx=$sizex/" ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini
    sed -i "s/^ny=[0-9]*$/ny=$size/" ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini
    sed -i "s/^nz=[0-9]*$/nz=$size/" ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini

    echo "$size"
    cat ${BASE_DIR}/envs/adastra/setup.ini | grep nx
    cat ${BASE_DIR}/envs/adastra/setup.ini | grep ny
    cat ${BASE_DIR}/envs/adastra/setup.ini | grep nz
    #cat ${BASE_DIR}/envs/adastra/io_deisa.yml | grep size
    echo "------"
    sbatch --wait -o resultsSizeBench/Deisa/${SIMU_SIZE}/res$size.out ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/launcher.sh
    sbatch --wait -o resultsSizeBench/NoDeisa/${SIMU_SIZE}/res$size.out ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/launcher_noDeisa.sh
done



#cat envs/adastra/setup.ini
#cat envs/adastra/io_deisa.yml