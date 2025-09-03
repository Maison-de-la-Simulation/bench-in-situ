#!/bin/bash

#SBATCH --account=$ACTIVE_PROJECT
#SBATCH --output=$1Bench.out
#SBATCH --job-name=b-i-s_nd
#SBATCH --constraint=GENOA
##SBATCH --constraint=$1
#SBATCH --nodes=1
##SBATCH --exclusive
#SBATCH --time=06:00:00
##SBATCH --nodelist=c1155
##SBATCH --gpus-per-node=1

source ../modules.env
source ../modules$1.env

export MPICH_GPU_SUPPORT_ENABLED=1

BASE_DIR=${PWD}
WORKING_DIR=${BASE_DIR}/working_dir
SIMU_SIZE=16
CUBE_SIZE=64
LAUNCHER_FILE="template_launcher.sh"
RESULT_DIR="results_bench_$1_$(date +%Y-%m-%d-%H-%M)"
RESULT_FILE=bench_output_$1_$(date +%Y-%m-%d-%H-%M).txt
FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")
FORMATED_CUBE_SIZE=$(printf "%03d" "$CUBE_SIZE")

mkdir -p ${WORKING_DIR} && cd ${WORKING_DIR}
rm *.h5
rm *.xmf
cd ..

declare -A local_domain_decomposition=(
    ['x']=0
    ['y']=0
    ['z']=0
)

echo $1 >> ${BASE_DIR}/${RESULT_FILE}
cat ${PWD}/../../../lib/pdi/pdi/VERSION >> ${BASE_DIR}/${RESULT_FILE}

for  ((CUBE_SIZE=${SMALL_CUBE_SIZE}; CUBE_SIZE<=${SMALL_CUBE_SIZE}; CUBE_SIZE*=2)); do
    FORMATED_CUBE_SIZE=$(printf "%03d" "$CUBE_SIZE")

    for SIMU_SIZE in "${!problem_subdivisions_small[@]}"; do
	FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")
        dd_value=${problem_subdivisions_small[$SIMU_SIZE]}
        echo "Key: $SIMU_SIZE ; Formated simu size: $FORMATED_SIMU_SIZE"

        JOB_GENERATED_DIR=${BASE_DIR}/${RESULT_DIR}/${FORMATED_CUBE_SIZE}/${FORMATED_SIMU_SIZE}
        mkdir -p ${JOB_GENERATED_DIR}
        GENERATED_SETUP_INI=${JOB_GENERATED_DIR}/setup.ini
        cp ${BASE_DIR}/../setup.ini ${GENERATED_SETUP_INI}
        echo $1 >> ${JOB_GENERATED_DIR}/metadata.dat
        cat ${PWD}/../../../lib/pdi/pdi/VERSION >> ${JOB_GENERATED_DIR}/metadata.dat
        GENERATED_LAUNCHER=${JOB_GENERATED_DIR}/launcher.sh
        cp ${LAUNCHER_FILE} ${GENERATED_LAUNCHER}
        cp ${BASE_DIR}/io_chkpt.yml ${JOB_GENERATED_DIR}/io_chkpt.yml
        cp ${BASE_DIR}/../modules.env ${JOB_GENERATED_DIR}/modules.env
        cp ${BASE_DIR}/../modules$1.env ${JOB_GENERATED_DIR}/modules$1.env

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
        let sizex=$CUBE_SIZE/${local_domain_decomposition['x']}
        let sizey=$CUBE_SIZE/${local_domain_decomposition['y']}
        let sizez=$CUBE_SIZE/${local_domain_decomposition['z']}
        echo "$sizex $sizey $sizez"
        
        # Change the discretisation between gpu
        sed -i "s/^nx=[0-9]*$/nx=$sizex/" ${GENERATED_SETUP_INI}
        sed -i "s/^ny=[0-9]*$/ny=$sizey/" ${GENERATED_SETUP_INI}
        sed -i "s/^nz=[0-9]*$/nz=$sizez/" ${GENERATED_SETUP_INI}

        sed -i "s/^mx=[0-9]*$/mx=${local_domain_decomposition['x']}/" ${GENERATED_SETUP_INI}
        sed -i "s/^my=[0-9]*$/my=${local_domain_decomposition['y']}/" ${GENERATED_SETUP_INI}
        sed -i "s/^mz=[0-9]*$/mz=${local_domain_decomposition['z']}/" ${GENERATED_SETUP_INI}

        sed -i "s/^#SBATCH --output=.*$/#SBATCH --output=res${FORMATED_SIMU_SIZE}N_%x_%j.out/" ${GENERATED_LAUNCHER}
        sed -i "s/^SIMU_SIZE=.*$/SIMU_SIZE=${SIMU_SIZE}/" ${GENERATED_LAUNCHER}
        sed -i "s/^SIM_PROC=.*$/SIM_PROC=${SIMU_SIZE}/" ${GENERATED_LAUNCHER}
        sed -i "s/^#SBATCH --account=.*$/#SBATCH --account=${ACTIVE_PROJECT}/" ${GENERATED_LAUNCHER}
        sed -i "s/^#SBATCH --constraint=.*$/#SBATCH --constraint=$1/" ${GENERATED_LAUNCHER}
        sed -i "s/modulesGPU.*$/modules$1.env/" ${GENERATED_LAUNCHER}
        sed -i "s|^\(JOB_GENERATED_DIR=\).*|\1${JOB_GENERATED_DIR}/|" ${GENERATED_LAUNCHER}
        if [ "$1" = "MI250" ]; then
            sed -i "s/^#SBATCH --cpus-per-task=.*$/#SBATCH --cpus-per-task=16/" ${GENERATED_LAUNCHER}
            if [ "$SIMU_SIZE" -gt 8 ]; then
                sed -i "s/^#SBATCH --nodes=.*$/#SBATCH --nodes=$((SIMU_SIZE / 8))/" ${GENERATED_LAUNCHER}
                sed -i "s/^SIM_NODES=.*$/SIM_NODES=$((SIMU_SIZE / 8))/" ${GENERATED_LAUNCHER}
            else
                sed -i "s/^#SBATCH --nodes=.*$/#SBATCH --nodes=1/" ${GENERATED_LAUNCHER}
                sed -i "s/^SIM_NODES=.*$/SIM_NODES=1/" ${GENERATED_LAUNCHER}
            fi
        elif [ "$1" = "MI300" ]; then
            sed -i "s/^#SBATCH --cpus-per-task=.*$/#SBATCH --cpus-per-task=24/" ${GENERATED_LAUNCHER}
            if [ "$SIMU_SIZE" -gt 4 ]; then
                sed -i "s/^#SBATCH --nodes=.*$/#SBATCH --nodes=$((SIMU_SIZE / 4))/" ${GENERATED_LAUNCHER}
                sed -i "s/^SIM_NODES=.*$/SIM_NODES=$((SIMU_SIZE / 4))/" ${GENERATED_LAUNCHER}
            else
                sed -i "s/^#SBATCH --nodes=.*$/#SBATCH --nodes=1/" ${GENERATED_LAUNCHER}
                sed -i "s/^SIM_NODES=.*$/SIM_NODES=1/" ${GENERATED_LAUNCHER}
            fi
        fi

        cat ${GENERATED_SETUP_INI} | grep nx
        cat ${GENERATED_SETUP_INI} | grep ny
        cat ${GENERATED_SETUP_INI} | grep nz

        echo "-----"
        cat ${GENERATED_SETUP_INI} | grep mx
        cat ${GENERATED_SETUP_INI} | grep my
        cat ${GENERATED_SETUP_INI} | grep mz

        echo "-----"
        cat ${GENERATED_LAUNCHER} | grep output=res
        cat ${GENERATED_LAUNCHER} | grep nodes=
        cat ${GENERATED_LAUNCHER} | grep ^SIMU_SIZE=
        cat ${GENERATED_LAUNCHER} | grep ^SIM_NODES=
        cat ${GENERATED_LAUNCHER} | grep SIM_PROC=
        cat ${GENERATED_LAUNCHER} | grep modulesMI
        cat ${GENERATED_LAUNCHER} | grep cpus-per-task
        cat ${GENERATED_LAUNCHER} | grep constraint
        cat ${GENERATED_LAUNCHER} | grep ^LOCAL_DIR=

        sbatch --wait -o ${JOB_GENERATED_DIR}/res${FORMATED_SIMU_SIZE}.out ${GENERATED_LAUNCHER}
        echo "----------------------------------------"
    done
done

mkdir -p ${BASE_DIR}/${RESULT_DIR} && cd ${BASE_DIR}/${RESULT_DIR} && grep -rw . -e "RESULT" >> ${BASE_DIR}/${RESULT_FILE}
rm -rf ${WORKING_DIR}
