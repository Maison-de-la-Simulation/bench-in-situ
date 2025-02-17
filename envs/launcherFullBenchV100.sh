#!/bin/bash

########################################################################################
## EXAMPLE OF USE: BENCHMARK FOR STRONG SCALING FOR A GIVEN SIZE
##
##  1) To run the scripts of bench in the default directory "jeanzay_V100_16G"
##
##          bash launcherFullBench.sh -gpu=V100 -gpuMEM=16G
##
##  2) To run the scripts of the bench in the directory "jeanzay_V100_16G_CUBESIZE_256"
##
##          bash launcherFullBench.sh -gpu=V100 -gpuMEM=16G -cubesize=256
##
## Remark: - The bench script are not launched in the same time.
##         Each bench script wait that the previous one is finish
##         - The value of gpuMEM=16G for the moment. (need to add the 32G case for JZ)
#########################################################################################
for i in "$@"
do
case $i in
    -gpu=*)
    GPU_ARCH="${i#*=}"
    ;;
    -gpuMEM=*)
    GPU_MEM="${i#*=}"
    ;;
    -cubesize=*)
    CUBE_SIZE="${i#*=}"
    ;;
#    -simusize=*)
#    SIMU_SIZE="${i#*=}"
#    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
    # unknown option
    ;;
esac
done

declare -A tab_repart=(
    ['2']=" 2 1 1 "
    ['4']=" 2 2 1 "
    ['8']=" 2 2 2 "
#    ['16']=" 4 2 2 "
#    ['32']=" 4 4 2 "
#    ['64']=" 4 4 4 "
#    ['128']=" 4 4 8 "
)


echo GPU_ARCH = ${GPU_ARCH}
echo GPU_MEM = ${GPU_MEM}
echo CUBE_SIZE = ${CUBE_SIZE}
#echo SIMU_SIZE = ${SIMU_SIZE}

NODES_ARCH_SUPER_FRIEND="_"${GPU_ARCH}

if [ "${GPU_ARCH}" == "" ]; then
    echo "Error -gpu is not given."
    exit 1
fi
if [ "${GPU_ARCH}" == "V100" ]; then
    if [ "${GPU_MEM}" == "" ]; then
        echo "Error: The value -gpuMEM msut be 16G or 32G on JZ"
        exit 1
    fi
    if [ "${GPU_MEM}" != "16G" ]; then
        echo "Error: The value -gpuMEM msut be 16G for the moment"
        exit 1
    fi
fi

if [ "${CUBE_SIZE}" == "" ]; then
    PDI_MHD_NODES_ARCH=${NODES_ARCH_SUPER_FRIEND}"_"${GPU_MEM}
else
    PDI_MHD_NODES_ARCH=${NODES_ARCH_SUPER_FRIEND}"_"${GPU_MEM}"_CUBESIZE_"${CUBE_SIZE}
fi


## Add test for the directory existence
echo NODES_ARCH_SUPER_FRIEND = ${NODES_ARCH_SUPER_FRIEND}
echo PDI_MHD_NODES_ARCH = ${PDI_MHD_NODES_ARCH}
exit 1
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
BASE_DIR=${WORK}/numpex/bench-in-situ_pdi_1_8
WORKING_DIR=${BASE_DIR}/working_dir${NODES_ARCH_SUPER_FRIEND}

cd ${BASE_DIR}
mkdir -p resultsBench${PDI_MHD_NODES_ARCH}/NoDeisa/${CUBE_SIZE}

SIMU_SIZE=1
lastjobid=$(sbatch --parsable -o resultsBench${PDI_MHD_NODES_ARCH}/NoDeisa/${CUBE_SIZE}/res${SIMU_SIZE}_.out ${BASE_DIR}/envs/jeanzay${PDI_MHD_NODES_ARCH}/${SIMU_SIZE}/launcher_noDeisa.sh)
echo "Key: $SIMU_SIZE"

for SIMU_SIZE in "${!tab_repart[@]}"; do
        value=${tab_repart[$SIMU_SIZE]}
        echo "Key: $SIMU_SIZE"

        lastjobid=$(sbatch --parsable -d afterok:${lastjobid} -o resultsBench${PDI_MHD_NODES_ARCH}/NoDeisa/${CUBE_SIZE}/res${SIMU_SIZE}_.out ${BASE_DIR}/envs/jeanzay${PDI_MHD_NODES_ARCH}/${SIMU_SIZE}/launcher_noDeisa.sh)
done
