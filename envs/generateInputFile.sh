#!/bin/bash

for i in "$@"
do
case $i in
    -localDomain=*)
    ## If the cubesize corresponding to the size of the sub-domain, localDomain=1
    ## If the cubesize corresponding to the size of the global domain, localDomain=0
    LOCAL_DOMAIN="${i#*=}"
    ;;
    -inputdir=*)
    DIR_INPUT="${i#*=}"
    ;;
    -outputdir=*)
    DIR_OUTPUT="${i#*=}"
    ;;
    -date=*)
    DATE_SEND="${i#*=}"
    ;;
    -cubesize=*)
    CUBE_SIZE="${i#*=}"
    ;;
    -simusize=*)
    SIMU_SIZE="${i#*=}"
    ;;
    -launchdir=*)
    PDI_MHD_NODES_ARCH="${i#*=}"
    ;;
    -subdir=*)
    SUB_DIR="${i#*=}"
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
    # unknown option
    ;;
esac
done

# constant array
declare -A tab_repart=(
    ['1']=" 1 1 1 "
    ['2']=" 2 1 1 "
    ['4']=" 2 2 1 "
    ['8']=" 2 2 2 "
    ['16']=" 4 2 2 "
    ['32']=" 4 4 2 "
    ['64']=" 4 4 4 "
    ['128']=" 8 4 4 "
    ['256']=" 8 8 4 "
)

# Tableau associatif pour les valeurs x, y, z
declare -A tab_nxyz=(
    ['x']=0
    ['y']=0
    ['z']=0
)


echo "LOCAL_DOMAIN=${LOCAL_DOMAIN}"
if [ "${LOCAL_DOMAIN}" == "" ]; then
    echo "cube size corresponds to the global size"
    LOCAL_DOMAIN=1
elif [ "${LOCAL_DOMAIN}" == "1" ]; then
    echo "cube size corresponds to the global size"
elif [ "${LOCAL_DOMAIN}" == "0" ]; then
    echo "cube size corresponds to the size of the sub domains."
else
    echo "The -localDomain=1 or -localDomain=0"
    echo "If localDomain is not given, we suppose -localDomain=false by default."
    exit 1
fi

if [ "${SIMU_SIZE}" == "" ]; then
    if [ "${CUBE_SIZE}" != "" ]; then
        echo "error no simu_size are given"
        exit 1
    fi
fi

if [ "${CUBE_SIZE}" == "" ]; then
    if [ "${SIMU_SIZE}" != "" ]; then
        echo "error no cube_size are given"
        exit 1
    fi
fi

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

echo "tab_nxyz: x=${tab_nxyz['x']} y=${tab_nxyz['y']} z=${tab_nxyz['z']}"
if [ "${LOCAL_DOMAIN}" == "1" ]; then
    let sizex=$CUBE_SIZE
    let sizey=$CUBE_SIZE
    let sizez=$CUBE_SIZE
else
    let sizex=$CUBE_SIZE/${tab_nxyz['x']}
    let sizey=$CUBE_SIZE/${tab_nxyz['y']}
    let sizez=$CUBE_SIZE/${tab_nxyz['z']}
fi
echo "$sizex $sizey $sizez"

if [ ! -d "${DIR_INPUT}/${SIMU_SIZE}" ]; then
    echo "The directory ${DIR_INPUT}/${SIMU_SIZE} doesn't exist."
    exit 1
fi

FILE_SETUPINI=${DIR_INPUT}/${SIMU_SIZE}/setup.ini
if [ ! -f "${FILE_SETUPINI}" ]; then
    echo "File ${FILE_SETUPINI} doesn't exist"
    exit 1
fi

FILE_LAUNCHER_NODEISA=${DIR_INPUT}/${SIMU_SIZE}/launcher_noDeisa.sh
if [ ! -f "${FILE_LAUNCHER_NODEISA}" ]; then
    echo "File ${FILE_LAUNCHER_NODEISA} doesn't exist"
    exit 1
fi

cd ${DIR_OUTPUT}
mkdir ${SUB_DIR}
cd ${SUB_DIR}

cp ../../${FILE_SETUPINI} .
cp ../../${FILE_LAUNCHER_NODEISA} .

sed_param=s/^PDI_MHD_NODES_ARCH=.*/PDI_MHD_NODES_ARCH=\"$PDI_MHD_NODES_ARCH\"/g
echo ${sed_param}
sed -i "$sed_param" launcher_noDeisa.sh

sed -i "s/^nx=[0-9]*$/nx=$sizex/" setup.ini
sed -i "s/^ny=[0-9]*$/ny=$sizey/" setup.ini
sed -i "s/^nz=[0-9]*$/nz=$sizez/" setup.ini

sed -i "s/^mx=[0-9]*$/mx=${tab_nxyz['x']}/" setup.ini
sed -i "s/^my=[0-9]*$/my=${tab_nxyz['y']}/" setup.ini
sed -i "s/^mz=[0-9]*$/mz=${tab_nxyz['z']}/" setup.ini
