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
SIMU_SIZE=1
CUBE_SIZE=64
LAUNCHER_FILE="launcher_noDeisa.sh"
RESULT_DIR="results_bench_$1_$(date +%Y-%m-%d-%H-%M)"
RESULT_FILE=bench_output_$1_$(date +%Y-%m-%d-%H-%M).txt
FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")
FORMATED_CUBE_SIZE=$(printf "%03d" "$CUBE_SIZE")

mkdir -p ${WORKING_DIR} && cd ${WORKING_DIR}
rm *.h5
rm *.xmf
cd ..

source ${BASE_DIR}/../topology.dat

declare -A subdivisions_of_iteration=(
    ['x']=0
    ['y']=0
    ['z']=0
)

echo $1 >> ${BASE_DIR}/${RESULT_FILE}
cat ${PWD}/../../../lib/pdi/pdi/VERSION >> ${BASE_DIR}/${RESULT_FILE}

for  ((CUBE_SIZE=${MINIMUM_CUBE_SIZE}; CUBE_SIZE<=${MAXIMUM_CUBE_SIZE}; CUBE_SIZE*=2)); do
    FORMATED_CUBE_SIZE=$(printf "%03d" "$CUBE_SIZE")

    for SIMU_SIZE in "${!problem_subdivisions[@]}"; do
	FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")
        value=${problem_subdivisions[$SIMU_SIZE]}
        echo "Key: $SIMU_SIZE ; Formated simu size: $FORMATED_SIMU_SIZE"

        ITERATION_FOLDER=${BASE_DIR}/${RESULT_DIR}/NoDeisa/${FORMATED_CUBE_SIZE}/${FORMATED_SIMU_SIZE}
        mkdir -p ${ITERATION_FOLDER}
        WHICH_SETUP=${ITERATION_FOLDER}/setup.ini
        cp ${BASE_DIR}/../setup.ini ${WHICH_SETUP}
        cat ${PWD}/../../../lib/pdi/pdi/VERSION >> ${ITERATION_FOLDER}/metadata.dat
        WHICH_LAUNCHER=${ITERATION_FOLDER}/launcher.sh
        cp ${LAUNCHER_FILE} ${WHICH_LAUNCHER}
        cp ${BASE_DIR}/../topology.dat ${ITERATION_FOLDER}/topology.dat
        cp ${BASE_DIR}/../io_chkpt.yml ${ITERATION_FOLDER}/io_chkpt.yml
        cp ${BASE_DIR}/../modules.env ${ITERATION_FOLDER}/modules.env
        cp ${BASE_DIR}/../modules$1.env ${ITERATION_FOLDER}/modules$1.env

        i=0
        
        for number in $value; do
            case $i in
                0) subdivisions_of_iteration['x']=$number ;;
                1) subdivisions_of_iteration['y']=$number ;;
                2) subdivisions_of_iteration['z']=$number ;;
            esac
        ((i++))
        done
        
        echo "subdivisions_of_iteration: x=${subdivisions_of_iteration['x']} y=${subdivisions_of_iteration['y']} z=${subdivisions_of_iteration['z']}"
        let sizex=$CUBE_SIZE/${subdivisions_of_iteration['x']}
        let sizey=$CUBE_SIZE/${subdivisions_of_iteration['y']}
        let sizez=$CUBE_SIZE/${subdivisions_of_iteration['z']}
        echo "$sizex $sizey $sizez"
        sed -i "s/^nx=[0-9]*$/nx=$sizex/" ${WHICH_SETUP}
        sed -i "s/^ny=[0-9]*$/ny=$sizey/" ${WHICH_SETUP}
        sed -i "s/^nz=[0-9]*$/nz=$sizez/" ${WHICH_SETUP}

        sed -i "s/^mx=[0-9]*$/mx=${subdivisions_of_iteration['x']}/" ${WHICH_SETUP}
        sed -i "s/^my=[0-9]*$/my=${subdivisions_of_iteration['y']}/" ${WHICH_SETUP}
        sed -i "s/^mz=[0-9]*$/mz=${subdivisions_of_iteration['z']}/" ${WHICH_SETUP}

        sed -i "s/^#SBATCH --output=res.*$/#SBATCH --output=res${SIMU_SIZE}N_%x_%j.out/" ${WHICH_LAUNCHER}
        sed -i "s/^SIMU_SIZE=.*$/SIMU_SIZE=${SIMU_SIZE}/" ${WHICH_LAUNCHER}
        sed -i "s/^SIM_PROC=.*$/SIM_PROC=${SIMU_SIZE}/" ${WHICH_LAUNCHER}
        sed -i "s/^#SBATCH --account=.*$/#SBATCH --account=${ACTIVE_PROJECT}/" ${WHICH_LAUNCHER}
        sed -i "s/^#SBATCH --constraint=.*$/#SBATCH --constraint=$1/" ${WHICH_LAUNCHER}
        sed -i "s/modulesMI.*$/modules$1_${FORMATED_SIMU_SIZE}.env/" ${WHICH_LAUNCHER}
        sed -i "s|^\(LOCAL_DIR=\).*|\1${ITERATION_FOLDER}/|" ${WHICH_LAUNCHER}
        if [ "$1" = "MI250" ]; then
            sed -i "s/^#SBATCH --cpus-per-task=.*$/#SBATCH --cpus-per-task=16/" ${WHICH_LAUNCHER}
            if [ "$SIMU_SIZE" -gt 8 ]; then
                sed -i "s/^#SBATCH --nodes=.*$/#SBATCH --nodes=$((SIMU_SIZE / 8))/" ${WHICH_LAUNCHER}
                sed -i "s/^SIM_NODES=.*$/SIM_NODES=$((SIMU_SIZE / 8))/" ${WHICH_LAUNCHER}
            else
                sed -i "s/^#SBATCH --nodes=.*$/#SBATCH --nodes=1/" ${WHICH_LAUNCHER}
                sed -i "s/^SIM_NODES=.*$/SIM_NODES=1/" ${WHICH_LAUNCHER}
            fi
        elif [ "$1" = "MI300" ]; then
            sed -i "s/^#SBATCH --cpus-per-task=.*$/#SBATCH --cpus-per-task=24/" ${WHICH_LAUNCHER}
            if [ "$SIMU_SIZE" -gt 4 ]; then
                sed -i "s/^#SBATCH --nodes=.*$/#SBATCH --nodes=$((SIMU_SIZE / 4))/" ${WHICH_LAUNCHER}
                sed -i "s/^SIM_NODES=.*$/SIM_NODES=$((SIMU_SIZE / 4))/" ${WHICH_LAUNCHER}
            else
                sed -i "s/^#SBATCH --nodes=.*$/#SBATCH --nodes=1/" ${WHICH_LAUNCHER}
                sed -i "s/^SIM_NODES=.*$/SIM_NODES=1/" ${WHICH_LAUNCHER}
            fi
        fi

        cat ${WHICH_SETUP} | grep nx
        cat ${WHICH_SETUP} | grep ny
        cat ${WHICH_SETUP} | grep nz

        echo "-----"
        cat ${WHICH_SETUP} | grep mx
        cat ${WHICH_SETUP} | grep my
        cat ${WHICH_SETUP} | grep mz

        echo "-----"
        cat ${WHICH_LAUNCHER} | grep output=res
        cat ${WHICH_LAUNCHER} | grep nodes=
        cat ${WHICH_LAUNCHER} | grep ^SIMU_SIZE=
        cat ${WHICH_LAUNCHER} | grep ^SIM_NODES=
        cat ${WHICH_LAUNCHER} | grep SIM_PROC=
        cat ${WHICH_LAUNCHER} | grep modulesMI
        cat ${WHICH_LAUNCHER} | grep cpus-per-task
        cat ${WHICH_LAUNCHER} | grep constraint
        cat ${WHICH_LAUNCHER} | grep ^LOCAL_DIR=

        sbatch --wait -o ${RESULT_DIR}/NoDeisa/${FORMATED_CUBE_SIZE}/res${FORMATED_SIMU_SIZE}.out ${WHICH_LAUNCHER}
        echo "----------------------------------------"
    done
done

mkdir -p ${BASE_DIR}/${RESULT_DIR}/NoDeisa && cd ${BASE_DIR}/${RESULT_DIR}/NoDeisa && grep -rw . -e "RESULT" >> ${BASE_DIR}/${RESULT_FILE}
rm -rf ${WORKING_DIR}
