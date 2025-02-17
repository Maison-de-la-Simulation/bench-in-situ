#!/bin/bash

########################################################################
## EXAMPLE OF USE: BENCHMARK FOR STRONG SCALING FOR A GIVEN SIZE
##
## bash launcherFullBench.sh -gpu=A100 -cubesize=256
##
## Remark: the bench script are not launched in the same time.
##         Each beanch script wait that the previous one is finish.
##
##         The directory jeanzay_A100_CUBESIZE_256 need to exist before.
#########################################################################
for i in "$@"
do
case $i in
    -gpu=*)
    GPU_ARCH="${i#*=}"
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
echo GPU_ARCH = ${GPU_ARCH}
echo CUBE_SIZE = ${CUBE_SIZE}

declare -A tab_repart=(
    ['2']=" 2 1 1 "
    ['4']=" 2 2 1 "
    ['8']=" 2 2 2 "
#    ['16']=" 4 2 2 "
#    ['32']=" 4 4 2 "
#    ['64']=" 4 4 4 "
#    ['128']=" 4 4 8 "
)


NODES_ARCH_SUPER_FRIEND="_"${GPU_ARCH}
PDI_MHD_NODES_ARCH=${NODES_ARCH_SUPER_FRIEND}"_CUBESIZE_"${CUBE_SIZE}

echo NODES_ARCH_SUPER_FRIEND = ${NODES_ARCH_SUPER_FRIEND}
echo PDI_MHD_NODES_ARCH = ${PDI_MHD_NODES_ARCH}

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
