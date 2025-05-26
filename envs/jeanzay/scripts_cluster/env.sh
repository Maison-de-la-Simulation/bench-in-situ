SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
### cluster information
CLUSTER_NAME=jeanzay
GPU_ARCH=H100

### root of "bench-in-situ" repository
BENCH_ROOT_DIR=${SCRIPT_DIR}/../../..

### WORKING_DIR
if [ "${GPU_ARCH}" != "" ]; then
	WORKING_DIR=${SCRIPT_DIR}/../working_dir_${GPU_ARCH}
else
	WORKING_DIR=${SCRIPT_DIR}/../working_dir
fi

### SIMULATION CODE
BUILD_SIMULATION_DIR=${WORKING_DIR}/sim/build
MAIN_EXE_DIR=${BUILD_SIMULATION_DIR}              ## executable main are in this directory

### PDI
#PDI_DIR=${BENCH_ROOT_DIR}/lib/pdi
PDI_BUILD_DIR=${WORKING_DIR}/pdi/build
PDI_INSTALL_DIR=${WORKING_DIR}/pdi/install

### python
PYTHON_ENV=deisa
PYTHON_WORKING_DIR="output/tmp"
PYTHON_VERSION=${PYTHON_VERSION:-"311"}
### dask
SCHEFILE=${WORKING_DIR}/scheduler.json
DASK_DISTRIBUTED__COMM__UCX__INFINIBAND="False"

mkdir -p ${WORKING_DIR}
mkdir -p ${WORKING_DIR}/${PYTHON_ENV}
mkdir -p ${PDI_BUILD_DIR}

print_env() {
	echo "============="
	echo "SCRIPT_DIR=${SCRIPT_DIR}"
	echo "WORKING_DIR=${WORKING_DIR}"
	echo "SCHEFILE=${SCHEFILE}"
	echo "BUILD_SIMULATION_DIR=${BUILD_SIMULATION_DIR}"
	echo "PYTHON_ENV=${PYTHON_ENV}"
	echo "PYTHON_VERSION=${PYTHON_VERSION}"
	echo "DASK_DISTRIBUTED__COMM__UCX__INFINIBAND=${DASK_DISTRIBUTED__COMM__UCX__INFINIBAND}"
	echo "============="
}

# if [ "${PYTHON_DEISA}" == "ON" ]; then
#     #expand aliases defined in the shell ~/.profile and ~/.bashrc
#     shopt -s expand_aliases
#     [ -f ~/.profile ] && source ~/.profile
#     [ -f ~/.bashrc ] && source ~/.bashrc
# fi
