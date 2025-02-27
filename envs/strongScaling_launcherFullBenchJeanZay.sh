#!/bin/bash

########################################################################################
## EXAMPLE OF USE: BENCHMARK FOR STRONG SCALING FOR A GIVEN SIZE
##
##  1) To run the scripts of bench in the default directory "jeanzay_V100_16G"
##
##          bash launcherFullBenchJZ.sh -case=strongScaling -gpu=V100 -gpuMEM=16G -cubesize=16 -fsimusize=file_nbgpu.txt
##
##  2) To run the scripts of the bench in the directory "jeanzay_A100
##
##          bash launcherFullBenchJZ.sh -case=strongScaling -gpu=A100 -cubesize=16 -fsimusize=file_nbgpu.txt
##
##  In the file_nbgpu.txt, each line correspond to the number of gpu that we want to test (see example file_nbgpu.txt)
##
##
##
## Remark: - The bench script are not launched in the same time.
##         Each bench script wait that the previous one is finish
##         - The value of gpuMEM=16G for the moment. (need to add the 32G case for JZ)
#########################################################################################

for i in "$@"
do
case $i in
    -case=*)
    CASE_BATCH="${i#*=}"
    ;;
    -gpu=*)
    GPU_ARCH="${i#*=}"
    ;;
    -gpuMEM=*)
    GPU_MEM="${i#*=}"
    ;;
    -cubesize=*)
    CUBE_SIZE="${i#*=}"
    ;;
    -simusize=*)
    SIMU_SIZE="${i#*=}"
    ;;
    ## read from file the parameter
    -fsimusize=*)
    FILE_SIMU_SIZE="${i#*=}"
    echo "FILE_SIMU_SIZE=${FILE_SIMU_SIZE}"
    mapfile -d '' array1 < "${FILE_SIMU_SIZE}"
    declare -p array1
    FILE_WITH_SIMU_SIZE="TRUE"
    SIMU_SIZEs=()
    for TMP_VAR in "${!array1[@]}"; do
        value=${array1[$TMP_VAR]}
        echo "value=$value"
        for number in $value; do
            echo "number=$number"
            SIMU_SIZEs+=( $number )
        done
    done
    ;;
    ## read from file the parameter
    -fcubesize=*)
    mapfile -d '' array2 < "${i#*=}"
    FILE_WITH_CUBE_SIZE_FOR_ONE_GPU="TRUE"
    CUBE_SIZEs=()
    for TMP_VAR in "${!array2[@]}"; do
        value=${array2[$TMP_VAR]}
        echo "value=$value"
        for number in $value; do
            echo "number=$number"
            CUBE_SIZEs+=( $number )
        done
    done
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
    # unknown option
    ;;
esac
done

CLUSTER_NAME=jeanzay ## TO BE USE WHEN WE WANT TO USE FOR OTHERS CLUSTERS

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## CHECK ENTRY

if [ "${CASE_BATCH}" == "strongScaling" ]; then
    echo "strong scaling experiment"
    echo CUBE_SIZE = ${CUBE_SIZE}
else
    echo "-case must be equal to strongScaling"
    exit 1
fi

NODES_ARCH_SUPER_FRIEND="_"${GPU_ARCH}

if [ "${GPU_ARCH}" == "" ]; then
    echo "Error -gpu is not given."
    exit 1
elif [ "${GPU_ARCH}" == "A100" ]; then
    echo "launch in A100"
elif [ "${GPU_ARCH}" == "V100" ]; then
    if [ "${GPU_MEM}" == "" ]; then
        echo "Error: The value -gpuMEM must be 16G or 32G on JZ"
        exit 1
    fi
    if [ "${GPU_MEM}" != "16G" ] && [ "${GPU_MEM}" != "32G" ]; then
        echo "Error: The value -gpuMEM must be 16G or 32G on JZ for the moment"
        exit 1
    fi
else
    echo "gpu arch must be V100 or A100"
    exit 1
fi


#################################################################################
# option a rajouter: BASE_DIR, NAME OF LAUNCHER, RESULT_DIR
# ecriture des resultats
# script envoyer pour l'instant dans envs (pas terrible)
#################################################################################
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
echo "SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
## need to be lunch in jeanzay the directory envs
BASE_DIR=${SCRIPT_DIR}/..
WORKING_DIR=${BASE_DIR}/working_dir${NODES_ARCH_SUPER_FRIEND}
DATE_SEND=$(date +"%Y%m%d_%H%M%S_%4N")


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## BEGIN CREATION OF THE DIRECTORY FOR THE BENCH

ADD_INFO_ARCH=""
if [ "${GPU_ARCH}" == "V100" ]; then
    ADD_INFO_ARCH="_"${GPU_MEM}
fi

if [ "${CUBE_SIZE}" == "" ]; then
    echo "Error: for a strong scaling experiment, we need to give the global size in each direction of the domain."
    echo "Give the global size with -cubesize"
    exit 1
else
    PDI_MHD_NODES_ARCH=${NODES_ARCH_SUPER_FRIEND}${ADD_INFO_ARCH}"_CUBESIZE_"${CUBE_SIZE}

    INPUT_DIR="jeanzay"${NODES_ARCH_SUPER_FRIEND}${ADD_INFO_ARCH}
    OUTPUT_DIR="jeanzay"${PDI_MHD_NODES_ARCH}
    ## creation of the directory
    mkdir ${SCRIPT_DIR}/${OUTPUT_DIR}
    echo "create file ${SCRIPT_DIR}/${OUTPUT_DIR}"

    cp ${INPUT_DIR}/* ${OUTPUT_DIR}
    ## change the directory of the snapshot
    sed_snapshot=s:bench-in-situ:bench-in-situ_$DATE_SEND:g
    sed -i "$sed_snapshot" ${OUTPUT_DIR}/modules.env

    declare -p SIMU_SIZEs
    for SIMU_SIZE in "${SIMU_SIZEs[@]}"; do
        echo " je lance le batch for ${SIMU_SIZE}"
        bash generateInputFile.sh -inputdir=${INPUT_DIR} -outputdir=${OUTPUT_DIR} \
            -cubesize=${CUBE_SIZE} -simusize=${SIMU_SIZE} -date=${DATA_SEND} \
            -launchdir=${PDI_MHD_NODES_ARCH} -subdir=${SIMU_SIZE}
    done
fi

## END CREATION OF THE DIRECTORY FOR THE BENCH
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## BEGIN LAUNCH SCRIPT

## Add test for testing if the directory exists
echo "NODES_ARCH_SUPER_FRIEND=${NODES_ARCH_SUPER_FRIEND}"

## lunch the script
RESULT_DIR=resultsBench${DATE_SEND}${PDI_MHD_NODES_ARCH}/NoDeisa/${CUBE_SIZE}

echo "RESULT_DIR=${RESULT_DIR}"
sleep 1 ## TO HAVE DATE BE CORRECT IN SECOND

cd ${BASE_DIR}
mkdir -p ${RESULT_DIR}

for ii in "${!SIMU_SIZEs[@]}"; do

    SIMU_SIZE=${SIMU_SIZEs[$ii]};
    echo "create batch =${ii} for SIMU_SIZE=${SIMU_SIZE}"

    if [ $ii == 0 ]; then
        lastjobid=$(sbatch --parsable -o ${RESULT_DIR}/res${SIMU_SIZE}_.out ${BASE_DIR}/envs/jeanzay${PDI_MHD_NODES_ARCH}/${SIMU_SIZE}/launcher_noDeisa.sh)
    else
        lastjobid=$(sbatch --parsable -d afterok:${lastjobid} -o ${RESULT_DIR}/res${SIMU_SIZE}_.out ${BASE_DIR}/envs/jeanzay${PDI_MHD_NODES_ARCH}/${SIMU_SIZE}/launcher_noDeisa.sh)
    fi

done

## END LAUNCH SCRIPT
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
