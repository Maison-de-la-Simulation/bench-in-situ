#!/bin/bash
for i in "$@"
do
case $i in
    -gpu=*)
    GPU_ARCH="${i#*=}"
    ;;
    -cluster=*)
    CLUSTER_NAME="${i#*=}"
    ;;
    -deisa=*)
    PYTHON_DEISA="${i#*=}"
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
    # unknown option
    ;;
esac
done

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
cd ${SCRIPT_DIR}
source env.sh -gpu=${GPU_ARCH} -cluster=${CLUSTER_NAME}

## if cluster add the modules on the compilation
# if [ "${GPU_ARCH}" != "" ]; then
#     if [ "${GPU_ARCH}" == "V100" ]; then
#         source ${SCRIPT_DIR}/../envs/jeanzay/modules_V100.env
#         module list
#     elif [ "${GPU_ARCH}" == "A100" ]; then
#         source ${SCRIPT_DIR}/../envs/jeanzay/modules_A100.env
#         module list
#     elif [ "${GPU_ARCH}" == "H100" ]; then
#         source ${SCRIPT_DIR}/../envs/jeanzay/modules_H100.env
#         module list
#     fi
# fi

if [ "${PYTHON_DEISA}" == "ON" ];then
    ./offline_python_env_setup.sh
fi
./build_pdi.sh -gpu=${GPU_ARCH} -cluster=${CLUSTER_NAME}
./build_simulation.sh -gpu=${GPU_ARCH} -cluster=${CLUSTER_NAME}

# copy files to working directory
cd ..
if [ "${PYTHON_DEISA}" == "ON" ];then
    cp deisa_deps_py* ${WORKING_DIR}
fi
cp -r in-situ ${WORKING_DIR}
#cp io.yml io_chkpt.yml io_deisa.yml ${WORKING_DIR}

