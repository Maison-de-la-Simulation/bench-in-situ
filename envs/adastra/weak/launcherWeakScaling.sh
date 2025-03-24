#!/bin/bash
#SBATCH --account=$ACTIVE_PROJECT
#sbatch --output=$1WeakScaling.out
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
WHICH_LAUNCHER="launcher_noDeisa.sh"
RESULT_DIR="results_weakScaling_$1_$(date +%Y-%m-%d-%H-%M)"
RESULT_FILE=bench_weakScaling_$1_$(date +%Y-%m-%d-%H-%M).txt
FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")
FORMATED_CUBE_SIZE=$(printf "%03d" "$CUBE_SIZE")

mkdir -p ${WORKING_DIR} && cd ${WORKING_DIR}
rm *.h5
rm *.xmf
cd ..

declare -A problem_subdivisions=(
    ['1']=" 1 1 1 "
    ['2']=" 2 1 1 "
    ['4']=" 2 2 1 "
    ['8']=" 2 2 2 "
    ['16']=" 4 2 2 "
)

declare -A subdivisions_of_iteration=(
    ['x']=0
    ['y']=0
    ['z']=0
)

grep -v "##*" -rw ${BASE_DIR}/${WHICH_LAUNCHER} | grep -e "#SBATCH --constraint=" >> ${BASE_DIR}/${RESULT_FILE}
cat ${PWD}/../../../lib/pdi/pdi/VERSION >> ${BASE_DIR}/${RESULT_FILE}

for SIMU_SIZE in "${!problem_subdivisions[@]}"; do
    FORMATED_SIMU_SIZE=$(printf "%03d" "$SIMU_SIZE")
    CUBE_SIZE=512
    value=${problem_subdivisions[$SIMU_SIZE]}
    echo "Key: $SIMU_SIZE, CUBE_SIZE: $CUBE_SIZE"

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
    echo "512 512 512"
    sed -i "s/^nx=[0-9]*$/nx=512/" ${BASE_DIR}/../setup.ini
    sed -i "s/^ny=[0-9]*$/ny=512/" ${BASE_DIR}/../setup.ini
    sed -i "s/^nz=[0-9]*$/nz=512/" ${BASE_DIR}/../setup.ini

    sed -i "s/^mx=[0-9]*$/mx=${subdivisions_of_iteration['x']}/" ${BASE_DIR}/../setup.ini
    sed -i "s/^my=[0-9]*$/my=${subdivisions_of_iteration['y']}/" ${BASE_DIR}/../setup.ini
    sed -i "s/^mz=[0-9]*$/mz=${subdivisions_of_iteration['z']}/" ${BASE_DIR}/../setup.ini

    sed -i "s/^#SBATCH --output=res.*$/#SBATCH --output=res${SIMU_SIZE}N_%x_%j.out/" ${BASE_DIR}/${WHICH_LAUNCHER}
    sed -i "s/^SIMU_SIZE=.*$/SIMU_SIZE=${SIMU_SIZE}/" ${BASE_DIR}/${WHICH_LAUNCHER}
    sed -i "s/^SIM_PROC=.*$/SIM_PROC=${SIMU_SIZE}/" ${BASE_DIR}/${WHICH_LAUNCHER}
    sed -i "s/^#SBATCH --account=.*$/#SBATCH --account=${ACTIVE_PROJECT}/" ${BASE_DIR}/${WHICH_LAUNCHER}
    sed -i "s/^#SBATCH --constraint=.*$/#SBATCH --constraint=$1/" ${BASE_DIR}/${WHICH_LAUNCHER}
    sed -i "s/modulesMI.*$/modules$1.env/" ${BASE_DIR}/${WHICH_LAUNCHER}
    if [ "$1" = "MI250" ]; then
        sed -i "s/^#SBATCH --cpus-per-task=.*$/#SBATCH --cpus-per-task=16/" ${BASE_DIR}/${WHICH_LAUNCHER}
        sed -i "s/^#SBATCH --nodes=.*$/#SBATCH --nodes=$((SIMU_SIZE / 8))/" ${BASE_DIR}/${WHICH_LAUNCHER}
        sed -i "s/^SIM_NODES=.*$/SIM_NODES=$((SIMU_SIZE / 8))/" ${BASE_DIR}/${WHICH_LAUNCHER}
    elif [ "$1" = "MI300" ]; then
        sed -i "s/^#SBATCH --cpus-per-task=.*$/#SBATCH --cpus-per-task=24/" ${BASE_DIR}/${WHICH_LAUNCHER}
        sed -i "s/^#SBATCH --nodes=.*$/#SBATCH --nodes=$((SIMU_SIZE / 4))/" ${BASE_DIR}/${WHICH_LAUNCHER}
        sed -i "s/^SIM_NODES=.*$/SIM_NODES=$((SIMU_SIZE / 4))/" ${BASE_DIR}/${WHICH_LAUNCHER}
    fi

    cat ${BASE_DIR}/../setup.ini | grep nx
    cat ${BASE_DIR}/../setup.ini | grep ny
    cat ${BASE_DIR}/../setup.ini | grep nz

    echo "-----"
    cat ${BASE_DIR}/../setup.ini | grep mx
    cat ${BASE_DIR}/../setup.ini | grep my
    cat ${BASE_DIR}/../setup.ini | grep mz

    echo "-----"
    cat ${BASE_DIR}/${WHICH_LAUNCHER} | grep output=res
    cat ${BASE_DIR}/${WHICH_LAUNCHER} | grep nodes=
    cat ${BASE_DIR}/${WHICH_LAUNCHER} | grep ^SIMU_SIZE=
    cat ${BASE_DIR}/${WHICH_LAUNCHER} | grep ^SIM_NODES=
    cat ${BASE_DIR}/${WHICH_LAUNCHER} | grep SIM_PROC=
    cat ${BASE_DIR}/${WHICH_LAUNCHER} | grep modulesMI
    cat ${BASE_DIR}/${WHICH_LAUNCHER} | grep cpus-per-task
    cat ${BASE_DIR}/${WHICH_LAUNCHER} | grep constraint

    sbatch --wait -o ${RESULT_DIR}/NoDeisa/${FORMATED_CUBE_SIZE}/res${FORMATED_SIMU_SIZE}.out ${BASE_DIR}/${WHICH_LAUNCHER}
    echo "----------------------------------------"
done

mkdir -p ${BASE_DIR}/${RESULT_DIR}/NoDeisa && cd ${BASE_DIR}/${RESULT_DIR}/NoDeisa && grep -rw . -e "RESULT" >> ${BASE_DIR}/${RESULT_FILE}
rm -rf ${WORKING_DIR}
