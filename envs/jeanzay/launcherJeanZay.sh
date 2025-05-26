#!/bin/bash

# Get the parameters of the launcher
source param_launcher.ini


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
echo "SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

## Needs to be launched from envs/jeanzay
BENCH_ROOT_DIR=${SCRIPT_DIR}/../..
INIT_DIR=${SCRIPT_DIR}
WORKING_DIR=${INIT_DIR}/working_dir_${GPU_ARCH}
PDI_INSTALL_DIR=${WORKING_DIR}/pdi/install
MAIN_EXE_DIR=${WORKING_DIR}/sim/build
DATE_SEND=$(date +"%Y%m%d_%H%M%S")
YAML_FILE=io_chkpt.yml

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## BEGIN: CHECK INPUT PARAMETERS

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

## END: CHECK INPUT PARAMETERS
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## BEGIN: CREATING FOLDERS/FILES FOR LAUNCHERS GENERATED

ADD_INFO_ARCH=""
if [ "${GPU_ARCH}" == "V100" ]; then
    ADD_INFO_ARCH="_"${GPU_MEM}
fi

if [ "${CUBE_SIZE}" == "" ]; then
    echo "Error: for a scaling experiment, we need to give the global size in each direction of the domain."
    echo "Give the global size with -cubesize"
    exit 1
fi

SIMU_ARCH=${GPU_ARCH}${ADD_INFO_ARCH}
## Launch the script
if [ "${CASE_BATCH}" == "strong" ]; then
    RESULT_DIR_NAME="strong/jeanzay_"${DATE_SEND}"_"${SIMU_ARCH}
elif [ "${CASE_BATCH}" == "weak" ]; then
    RESULT_DIR_NAME="weak/jeanzay_"${DATE_SEND}"_"${SIMU_ARCH}
else
    echo "error -case must be strong or weak"
    exit 1
fi

echo "RESULT_DIR_NAME=$RESULT_DIR_NAME"
FORMATED_CUBE_SIZE=$(printf "%03d" "$CUBE_SIZE")
RESULT_DIR=${INIT_DIR}/${RESULT_DIR_NAME}/${FORMATED_CUBE_SIZE}

## Creation of the directory
mkdir -p ${RESULT_DIR}
echo "RESULT_DIR=${RESULT_DIR}"

## Local domain decomposition
declare -A local_domain_decomposition=(
    ['x']=0
    ['y']=0
    ['z']=0
)

declare -p SIMU_SIZE_ARRAY
for SIMU_SIZE in "${SIMU_SIZE_ARRAY[@]}"; do
    FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")
    echo " creation of the files for the launcher script ${SIMU_SIZE}"

    # Set directory name of the launcher generated
    JOB_GENERATED_DIR=${RESULT_DIR}/nb_gpu_${FORMATED_SIMU_SIZE}

    mkdir -p ${JOB_GENERATED_DIR}
    echo "create file JOB_GENERATED_DIR=${JOB_GENERATED_DIR}"

    # Get domain decomposition(dd) for SIMU_SIZE
    dd_value=${domain_decomposition_for_each_simulation_size[$SIMU_SIZE]}
    # Index to get the domain decomposition of an axis
    axis_index=0

    # Get domain decomposition for each axis (x,y,z)
    for dd_axis in $dd_value; do
        case $axis_index in
            0) local_domain_decomposition['x']=$dd_axis ;;
            1) local_domain_decomposition['y']=$dd_axis ;;
            2) local_domain_decomposition['z']=$dd_axis ;;
        esac
        ((axis_index++))
    done

    echo "local_domain_decomposition: x=${local_domain_decomposition['x']} y=${local_domain_decomposition['y']} z=${local_domain_decomposition['z']}"
    if [ "${CASE_BATCH}" == "weak" ]; then
        let sizex=$CUBE_SIZE
        let sizey=$CUBE_SIZE
        let sizez=$CUBE_SIZE
    elif [ "${CASE_BATCH}" == "strong" ]; then
        let sizex=$CUBE_SIZE/${local_domain_decomposition['x']}
        let sizey=$CUBE_SIZE/${local_domain_decomposition['y']}
        let sizez=$CUBE_SIZE/${local_domain_decomposition['z']}
    fi

    ## Copy all files needed by the jobs
    cp ${INIT_DIR}/setup.ini               ${JOB_GENERATED_DIR}/setup.ini
    cp ${INIT_DIR}/modules_${GPU_ARCH}.env ${JOB_GENERATED_DIR}/modules_${GPU_ARCH}.env
    cp ${INIT_DIR}/${YAML_FILE}            ${JOB_GENERATED_DIR}/${YAML_FILE}
    cp ${INIT_DIR}/${LAUNCH_NAME}          ${JOB_GENERATED_DIR}/job.sh

    ## Update local setup.ini
    ## Set local filename
    GENERATED_SETUP_INI=${JOB_GENERATED_DIR}/setup.ini
    echo "JOB_SETUP_INI=$GENERATED_SETUP_INI"

    ## Change the discretisation between gpu
    sed -i "s/^nx=[0-9]*$/nx=$sizex/" ${GENERATED_SETUP_INI}
    sed -i "s/^ny=[0-9]*$/ny=$sizey/" ${GENERATED_SETUP_INI}
    sed -i "s/^nz=[0-9]*$/nz=$sizez/" ${GENERATED_SETUP_INI}

    sed -i "s/^mx=[0-9]*$/mx=${local_domain_decomposition['x']}/" ${GENERATED_SETUP_INI}
    sed -i "s/^my=[0-9]*$/my=${local_domain_decomposition['y']}/" ${GENERATED_SETUP_INI}
    sed -i "s/^mz=[0-9]*$/mz=${local_domain_decomposition['z']}/" ${GENERATED_SETUP_INI}

    ## Update local modules_${GPU_ARCH}.env
    ## Change the directory of the snapshot
    sed_snapshot=s:bench-in-situ:bench-in-situ_$DATE_SEND:g
    sed -i "$sed_snapshot" ${JOB_GENERATED_DIR}/modules_${GPU_ARCH}.env


    ## Update local launcher
    ## Set local filename
    GENERATED_LAUNCHER=${JOB_GENERATED_DIR}"/job.sh"
    echo "GENERATED_LAUNCHER=$GENERATED_LAUNCHER"

    ## Set global variables
    sed -i "s:^GPU_ARCH=.*$:GPU_ARCH=${GPU_ARCH}:g" ${GENERATED_LAUNCHER}
    sed -i "s:^YAML_FILE=.*$:YAML_FILE=$YAML_FILE:g" ${GENERATED_LAUNCHER}

    sed -i "s:%%PDI_INSTALL_DIR%%:${PDI_INSTALL_DIR}:g" ${GENERATED_LAUNCHER}
    sed -i "s:%%MAIN_EXE_DIR%%:${MAIN_EXE_DIR}:g" ${GENERATED_LAUNCHER}

    ## Change the output file name
    sed -i "s/^#SBATCH --output=res.*$/#SBATCH --output=res${SIMU_SIZE}_%x_%j.out/g" ${GENERATED_LAUNCHER}

    ## REMARK: THE VERSION WITH PARTITION gpu_p2 (V100 octo GPU) IS CURRENTLY MISSING

    ## Set account and partition constraint
    if [ "${GPU_ARCH}" == "V100" ]; then

        sed -i "s/^#SBATCH -A.*$/#SBATCH -A ${IDRPROJ}@v100/g" ${GENERATED_LAUNCHER}
        if [ "${GPU_MEM}" == "16G" ]; then
            sed -i "s/^#SBATCH -C.*$/#SBATCH -C v100\-16g/g" ${GENERATED_LAUNCHER}
        elif [ "${GPU_MEM}" == "32G" ]; then
            sed -i "s/^#SBATCH -C.*$/#SBATCH -C v100\-32g/g" ${GENERATED_LAUNCHER}
        fi
        sed -i "s|^#SBATCH --cpus-per-task=.*$|#SBATCH --cpus-per-task=10|g" ${GENERATED_LAUNCHER}

    elif [ "${GPU_ARCH}" = "A100" ]; then

        sed -i "s/^#SBATCH -A.*$/#SBATCH -A ${IDRPROJ}@a100/g" ${GENERATED_LAUNCHER}
        sed -i "s/^#SBATCH -C.*$/#SBATCH -C a100/g" ${GENERATED_LAUNCHER}
        sed -i "s|^#SBATCH --cpus-per-task=.*$|#SBATCH --cpus-per-task=8|g" ${GENERATED_LAUNCHER}

    elif [ "${GPU_ARCH}" = "H100" ]; then

        sed -i "s/^#SBATCH -A.*$/#SBATCH -A ${IDRPROJ}@h100/g" ${GENERATED_LAUNCHER}
        sed -i "s/^#SBATCH -C.*$/#SBATCH -C h100/g" ${GENERATED_LAUNCHER}
        sed -i "s|^#SBATCH --cpus-per-task=.*$|#SBATCH --cpus-per-task=24|g" ${GENERATED_LAUNCHER}

    fi

    ## Set nodes configurations
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
        sed -i "s/^##SBATCH --ntasks=.*$/##SBATCH --ntasks=$SIMU_SIZE/g" ${GENERATED_LAUNCHER}
        sed -i "s/^##SBATCH --nodes=.*$/#SBATCH --nodes=1/g" ${GENERATED_LAUNCHER}
        sed -i "s/^#SBATCH --ntasks-per-node=.*$/#SBATCH --ntasks-per-node=$SIMU_SIZE/g" ${GENERATED_LAUNCHER}
        sed -i "s/^#SBATCH --gres=gpu:.*$/#SBATCH --gres=gpu:$SIMU_SIZE/g" ${GENERATED_LAUNCHER}

        sed -i "s/^SIMU_SIZE=.*$/SIMU_SIZE=$SIMU_SIZE/g" ${GENERATED_LAUNCHER}
    else
        ## CASE MULTI NODE
        SIMU_NODES=$((SIMU_SIZE / NB_GPU_PER_NODES))

        sed -i "s/^##SBATCH --ntasks=.*$/#SBATCH --ntasks=$SIMU_SIZE/g" ${GENERATED_LAUNCHER}
        sed -i "s/^##SBATCH --nodes=.*$/#SBATCH --nodes=$SIMU_NODES/g" ${GENERATED_LAUNCHER}
        sed -i "s/^#SBATCH --ntasks-per-node=.*$/#SBATCH --ntasks-per-node=$NB_GPU_PER_NODES/g" ${GENERATED_LAUNCHER}
        sed -i "s/^#SBATCH --gres=gpu:.*$/#SBATCH --gres=gpu:$NB_GPU_PER_NODES/g" ${GENERATED_LAUNCHER}

        sed -i "s/^SIMU_SIZE=.*$/SIMU_SIZE=$SIMU_SIZE/g" ${GENERATED_LAUNCHER}
    fi
done

## END: CREATING FOLDERS/FILES FOR LAUNCHERS GENERATED
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## BEGIN: SUBMISSION OF JOBS GENERATED

## Add test for testing if the directory exists
echo "RESULT_DIR_NAME=${RESULT_DIR_NAME}"
echo "CUBE_SIZE=$CUBE_SIZE"
sleep 1 ## TO PREVENT DIFFERENT JOBS FROM WRITING IN THE SAME FOLDER

index_job=1
for SIMU_SIZE in "${SIMU_SIZE_ARRAY[@]}"; do
    FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")
    JOB_GENERATED_DIR=${RESULT_DIR}/nb_gpu_${FORMATED_SIMU_SIZE}


    cd ${JOB_GENERATED_DIR}

    echo "JOB_GENERATED_DIR=$JOB_GENERATED_DIR"
    if [ ! -d "$JOB_GENERATED_DIR" ]; then
        echo "$JOB_GENERATED_DIR does not exist."
        exit 1
    fi

    if [ $index_job -eq 1 ]; then
        lastjobid=$(sbatch --parsable -o res${FORMATED_SIMU_SIZE}.out job.sh)
    else
        lastjobid=$(sbatch --parsable -d afterok:${lastjobid} -o res${FORMATED_SIMU_SIZE}.out job.sh)
    fi

    echo "Create batch=${index_job} for SIMU_SIZE=${SIMU_SIZE}"
    ((index_job++))
done

## END: SUBMISSION OF JOBS GENERATED
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
