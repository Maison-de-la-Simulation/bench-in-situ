#!/bin/bash
#SBATCH --account=cad14985
#sbatch --output=ALLSizeBench.out
#SBATCH --job-name=b-i-s_nd
#SBATCH --constraint=GENOA
##SBATCH --constraint=MI250
#SBATCH --nodes=1
##SBATCH --exclusive
#SBATCH --time=20:00:00
##SBATCH --nodelist=c1155
##SBATCH --gpus-per-node=1

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
BASE_DIR=${PWD}
WORKING_DIR=${BASE_DIR}/working_dir
SIMU_SIZE=1
CUBE_SIZE=64
WHICH_LAUNCHER="launcher_noDeisa-SOY25_MPICH.sh"
RESULT_DIR=results_v4_ALL_MPICH
RESULT_FILE=ALL_MPICH-output_v4.txt

cd ${WORKING_DIR}
rm *.h5
rm *.xmf
cd ..

declare -A tab_repart=(
#    ['1']=" 1 1 1 "
#    ['2']=" 2 1 1 "
    ['4']=" 2 2 1 "
    ['8']=" 2 2 2 "
    ['16']=" 4 2 2 "
    ['32']=" 4 4 2 "
#    ['64']=" 4 4 4 "
#    ['128']=" 4 4 8 "
)

# Tableau associatif pour les valeurs x, y, z
declare -A tab_nxyz=(
    ['x']=0
    ['y']=0
    ['z']=0
)

grep -v "##*" -rw ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/${WHICH_LAUNCHER} | grep -e "#SBATCH --constraint=" >> ${BASE_DIR}/${RESULT_FILE}
cat ${PWD}/lib/pdi/pdi/VERSION >> ${BASE_DIR}/${RESULT_FILE}

for  ((CUBE_SIZE=64; CUBE_SIZE<=128; CUBE_SIZE+=64)); do
    # Boucle sur chaque clé du tableau associatif
    for SIMU_SIZE in "${!tab_repart[@]}"; do
        value=${tab_repart[$SIMU_SIZE]}
        echo "Key: $SIMU_SIZE"
        
        # Compteur pour assigner les valeurs aux clés x, y, z
        i=0
        
        # Lire les chiffres un par un
        for number in $value; do
            case $i in
                0) tab_nxyz['x']=$number ;;
                1) tab_nxyz['y']=$number ;;
                2) tab_nxyz['z']=$number ;;
            esac
        ((i++))
        done
        
        # Afficher les valeurs du tableau tab_nxyz
        echo "tab_nxyz: x=${tab_nxyz['x']} y=${tab_nxyz['y']} z=${tab_nxyz['z']}"
        let sizex=$CUBE_SIZE/${tab_nxyz['x']}
        let sizey=$CUBE_SIZE/${tab_nxyz['y']}
        let sizez=$CUBE_SIZE/${tab_nxyz['z']}
        echo "$sizex $sizey $sizez"
        sed -i "s/^nx=[0-9]*$/nx=$sizex/" ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini
        sed -i "s/^ny=[0-9]*$/ny=$sizey/" ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini
        sed -i "s/^nz=[0-9]*$/nz=$sizez/" ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini
        
        sed -i "s/^mx=[0-9]*$/mx=${tab_nxyz['x']}/" ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini
        sed -i "s/^my=[0-9]*$/my=${tab_nxyz['y']}/" ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini
        sed -i "s/^mz=[0-9]*$/mz=${tab_nxyz['z']}/" ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini

        cat ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini | grep nx
        cat ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini | grep ny
        cat ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini | grep nz

        echo "-----"
        cat ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini | grep mx
        cat ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini | grep my
        cat ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/setup.ini | grep mz
        sbatch --wait -o ${RESULT_DIR}/NoDeisa/${CUBE_SIZE}/res${SIMU_SIZE}.out ${BASE_DIR}/envs/adastra/${SIMU_SIZE}/${WHICH_LAUNCHER}
        echo "----------------------------------------"
    done
done

# mkdir -p ${BASE_DIR}/${RESULT_DIR}/NoDeisa && cd ${BASE_DIR}/${RESULT_DIR}/NoDeisa && grep -rw ${CUBE_SIZE} -e "RESULT" >> ${BASE_DIR}/ALL-${CUBE_SIZE}-output_v4.txt
mkdir -p ${BASE_DIR}/${RESULT_DIR}/NoDeisa && cd ${BASE_DIR}/${RESULT_DIR}/NoDeisa && grep -rw . -e "RESULT" >> ${BASE_DIR}/${RESULT_FILE}
