#!/bin/bash

# get the parameter of the launcher
source param_launcher.ini


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
echo "SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

## need to be launch in envs/jeanzay
BENCH_ROOT_DIR=${SCRIPT_DIR}/../..
INIT_DIR=${SCRIPT_DIR}
RESULT_ROOT_DIR=${INIT_DIR}
WORKING_DIR=${INIT_DIR}/working_dir_${GPU_ARCH}
PDI_INSTALL_DIR=${WORKING_DIR}/pdi/install
MAIN_EXE_DIR=${WORKING_DIR}/sim/build
DATE_SEND=$(date +"%Y%m%d_%H%M%S")
YAML_FILE=io_chkpt.yml

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## CHECK ENTRY

if [ "${CASE_BATCH}" == "strong" ]; then
    echo "strong scaling experiment"
    echo CUBE_SIZE = ${CUBE_SIZE}
elif [ "${CASE_BATCH}" == "weak" ]; then
    echo "weak scaling experiment"
    echo CUBE_SIZE = ${CUBE_SIZE}
else
    echo "-case must be equal to strong or weak"
    exit 1
fi

if [ "${GPU_ARCH}" == "" ]; then
    echo "Error -gpu is not given."
    exit 1
elif [ "${GPU_ARCH}" == "A100" ]; then
    echo "launch in A100"
elif [ "${GPU_ARCH}" == "H100" ]; then
    echo "launch in H100"
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

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## BEGIN CREATION OF THE DIRECTORY FOR THE BENCH

ADD_INFO_ARCH=""
if [ "${GPU_ARCH}" == "V100" ]; then
    ADD_INFO_ARCH="_"${GPU_MEM}
fi

if [ "${CUBE_SIZE}" == "" ]; then
    echo "Error: for a scaling experiment, we need to give the global size in each direction of the domain."
    echo "Give the global size with -cubesize"
    exit 1
fi

PDI_MHD_NODES_ARCH=${GPU_ARCH}${ADD_INFO_ARCH}
## lunch the script
if [ "${CASE_BATCH}" == "strong" ]; then
    RESULT_DIR_NAME="strong/jeanzay_"${DATE_SEND}"_"${PDI_MHD_NODES_ARCH}
elif [ "${CASE_BATCH}" == "weak" ]; then
    RESULT_DIR_NAME="weak/jeanzay_"${DATE_SEND}"_"${PDI_MHD_NODES_ARCH}
else
    echo "error -case must be strong or weak"
    exit 1
fi

echo "RESULT_DIR_NAME=$RESULT_DIR_NAME"
FORMATED_CUBE_SIZE=$(printf "%03d" "$CUBE_SIZE")
RESULT_DIR=${RESULT_ROOT_DIR}/${RESULT_DIR_NAME}/${FORMATED_CUBE_SIZE}

## creation of the directory
mkdir -p ${RESULT_DIR}
echo "create file RESULT_DIR=${RESULT_DIR}"

# Tableau associatif pour les valeurs x, y, z
declare -A tab_nxyz=(
    ['x']=0
    ['y']=0
    ['z']=0
)

declare -p SIMU_SIZEs
for SIMU_SIZE in "${SIMU_SIZEs[@]}"; do
    FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")
    echo " creation of input file for the launcher script ${SIMU_SIZE}"

    TEST_CASE_DIR=${RESULT_DIR}/nb_gpu_${FORMATED_SIMU_SIZE}

    mkdir -p ${TEST_CASE_DIR}
    echo "create file TEST_CASE_DIR=${TEST_CASE_DIR}"

    value=${tab_repart[$SIMU_SIZE]}
    # Compteur pour assigner les valeurs aux clés x, y, z
    index=0

    # Lire les chiffres un par un
    for number in $value; do
        case $index in
            0) tab_nxyz['x']=$number ;;
            1) tab_nxyz['y']=$number ;;
            2) tab_nxyz['z']=$number ;;
        esac
        ((index++))
    done

    echo "tab_nxyz: x=${tab_nxyz['x']} y=${tab_nxyz['y']} z=${tab_nxyz['z']}"
    if [ "${CASE_BATCH}" == "weak" ]; then
        let sizex=$CUBE_SIZE
        let sizey=$CUBE_SIZE
        let sizez=$CUBE_SIZE
    elif [ "${CASE_BATCH}" == "strong" ]; then
        let sizex=$CUBE_SIZE/${tab_nxyz['x']}
        let sizey=$CUBE_SIZE/${tab_nxyz['y']}
        let sizez=$CUBE_SIZE/${tab_nxyz['z']}
    fi

    cp ${INIT_DIR}/setup.ini ${TEST_CASE_DIR}/setup.ini
    cp ${INIT_DIR}/modules_${GPU_ARCH}.env ${TEST_CASE_DIR}/modules_${GPU_ARCH}.env
    cp ${INIT_DIR}/${YAML_FILE} ${TEST_CASE_DIR}/${YAML_FILE}
    cp ${INIT_DIR}/${LAUNCH_NAME} ${TEST_CASE_DIR}/job.sh

    ## change the discretisation between gpu
    sed -i "s/^nx=[0-9]*$/nx=$sizex/" ${TEST_CASE_DIR}/setup.ini
    sed -i "s/^ny=[0-9]*$/ny=$sizey/" ${TEST_CASE_DIR}/setup.ini
    sed -i "s/^nz=[0-9]*$/nz=$sizez/" ${TEST_CASE_DIR}/setup.ini

    sed -i "s/^mx=[0-9]*$/mx=${tab_nxyz['x']}/" ${TEST_CASE_DIR}/setup.ini
    sed -i "s/^my=[0-9]*$/my=${tab_nxyz['y']}/" ${TEST_CASE_DIR}/setup.ini
    sed -i "s/^mz=[0-9]*$/mz=${tab_nxyz['z']}/" ${TEST_CASE_DIR}/setup.ini

    ## change the directory of the snapshot
    sed_snapshot=s:bench-in-situ:bench-in-situ_$DATE_SEND:g
    sed -i "$sed_snapshot" ${TEST_CASE_DIR}/modules_${GPU_ARCH}.env

    LOCAL_LAUNCHER=${TEST_CASE_DIR}"/job.sh"
    echo "LOCAL_LAUNCHER=$LOCAL_LAUNCHER"
    ## Tranform local launcher

    ## set global variable
    sed -i "s:^GPU_ARCH=.*$:GPU_ARCH=${GPU_ARCH}:" ${LOCAL_LAUNCHER}
    sed -i "s:^YAML_FILE=.*$:YAML_FILE=$YAML_FILE:" ${LOCAL_LAUNCHER}

    ## change the output
    sed -i "s/^#SBATCH --output=res.*$/#SBATCH --output=res${SIMU_SIZE}_%x_%j.out/" ${LOCAL_LAUNCHER}

    ## REMARK: WE DON'T HAVE CODED THE VERSION WITH PARTITION gpu_p2 (V100 octo GPU)

    ## set account and constraint
    if [ "${GPU_ARCH}" == "V100" ]; then
        sed -i "s/^#SBATCH -A.*$/#SBATCH -A ${IDRPROJ}@v100/g" ${LOCAL_LAUNCHER}
        if [ "${GPU_MEM}" == "16G" ]; then
            sed -i "s/^#SBATCH -C.*$/#SBATCH -C v100\-16g/g" ${LOCAL_LAUNCHER}
        elif [ "${GPU_MEM}" == "32G" ]; then
            sed -i "s/^#SBATCH -C.*$/#SBATCH -C v100\-32g/g" ${LOCAL_LAUNCHER}
        fi
        sed -i "s|^#SBATCH --cpus-per-task=.*$|#SBATCH --cpus-per-task=10           # nombre de CPU par tache (1/4 du noeud 4-GPU V100 ici)|g" ${LOCAL_LAUNCHER}

    elif [ "${GPU_ARCH}" = "A100" ]; then

        sed -i "s/^#SBATCH -A.*$/#SBATCH -A ${IDRPROJ}@a100/g" ${LOCAL_LAUNCHER}
        sed -i "s/^#SBATCH -C.*$/#SBATCH -C a100/g" ${LOCAL_LAUNCHER}
        sed -i "s|^#SBATCH --cpus-per-task=.*$|#SBATCH --cpus-per-task=8           # nombre de CPU par tache pour gpu_p5 (1/8 du noeud 8-GPU A100)|g" ${LOCAL_LAUNCHER}

    elif [ "${GPU_ARCH}" = "H100" ]; then

        sed -i "s/^#SBATCH -A.*$/#SBATCH -A ${IDRPROJ}@h100/g" ${LOCAL_LAUNCHER}
        sed -i "s/^#SBATCH -C.*$/#SBATCH -C h100/g" ${LOCAL_LAUNCHER}
        sed -i "s|^#SBATCH --cpus-per-task=.*$|#SBATCH --cpus-per-task=24          # nombre de CPU par tache pour gpu_p6 (1/4 du noeud 4-GPU H100)|g" ${LOCAL_LAUNCHER}

    fi

    ## set nodes configurations
    ## CASE 1: PARTITION WITH 4 GPU PER NODES
    if [ "${GPU_ARCH}" == "V100" ] || [ "${GPU_ARCH}" == "H100" ]; then
        NB_GPU_PER_NODES=4
    ## CASE 2: PARTITION WITH 8 GPU PER NODES
    elif [ "${GPU_ARCH}" == "A100" ]; then
        NB_GPU_PER_NODES=8
    else
        echo "GPU_ARCH=$GPU_ARCH doesn't exist."
        exit 1
    fi

    if [ $((SIMU_SIZE % NB_GPU_PER_NODES)) != 0 ] && [ $SIMU_SIZE -gt $NB_GPU_PER_NODES ]; then
        echo "SIMU_SIZE=${SIMU_SIZE} must be a power of ${NB_GPU_PER_NODES} if SIMU_SIZE=${SIMU_SIZE} supperior to ${NB_GPU_PER_NODES}."
        exit 1
    fi

    if [ $SIMU_SIZE -le $NB_GPU_PER_NODES ]; then
        ## CASE ONE NODE
        sed -i "s/^##SBATCH --ntasks=.*$/##SBATCH --ntasks=$SIMU_SIZE/g" ${LOCAL_LAUNCHER}
        sed -i "s/^##SBATCH --nodes=.*$/#SBATCH --nodes=1/" ${LOCAL_LAUNCHER}
        sed -i "s/^#SBATCH --ntasks-per-node=.*$/#SBATCH --ntasks-per-node=$SIMU_SIZE/g" ${LOCAL_LAUNCHER}
        sed -i "s/^#SBATCH --gres=gpu:.*$/#SBATCH --gres=gpu:$SIMU_SIZE/g" ${LOCAL_LAUNCHER}

        sed -i "s/^SIMU_SIZE=.*$/SIMU_SIZE=$SIMU_SIZE/g" ${LOCAL_LAUNCHER}
    else
        ## CASE MULTI NODE
        SIMU_NODES=$((SIMU_SIZE / NB_GPU_PER_NODES))

        sed -i "s/^##SBATCH --ntasks=.*$/#SBATCH --ntasks=$SIMU_SIZE/g" ${LOCAL_LAUNCHER}
        sed -i "s/^##SBATCH --nodes=.*$/#SBATCH --nodes=$SIMU_NODES/g" ${LOCAL_LAUNCHER}
        sed -i "s/^#SBATCH --ntasks-per-node=.*$/#SBATCH --ntasks-per-node=$NB_GPU_PER_NODES/g" ${LOCAL_LAUNCHER}
        sed -i "s/^#SBATCH --gres=gpu:.*$/#SBATCH --gres=gpu:$NB_GPU_PER_NODES/g" ${LOCAL_LAUNCHER}

        sed -i "s/^SIMU_SIZE=.*$/SIMU_SIZE=$SIMU_SIZE/g" ${LOCAL_LAUNCHER}
    fi

    ## set global variable
done

## END CREATION OF THE DIRECTORY FOR THE BENCH
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## BEGIN LAUNCH SCRIPT

## Add test for testing if the directory exists
echo "RESULT_DIR=${RESULT_DIR}"
sleep 1 ## TO HAVE THE CORRECT DATE IN SECOND

for ii in "${!SIMU_SIZEs[@]}"; do

    SIMU_SIZE=${SIMU_SIZEs[$ii]};
    FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")
    LAUNCH_TEST_CASE_DIR=${RESULT_DIR}/nb_gpu_${FORMATED_SIMU_SIZE}

    echo "create batch=${ii} for SIMU_SIZE=${SIMU_SIZE}"

    cd ${LAUNCH_TEST_CASE_DIR}

    if [ $ii == 0 ]; then
        lastjobid=$(sbatch --parsable -o res${FORMATED_SIMU_SIZE}.out job.sh)
    else
        lastjobid=$(sbatch --parsable -d afterok:${lastjobid} -o res${FORMATED_SIMU_SIZE}.out job.sh)
    fi

done

## END LAUNCH SCRIPT
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
