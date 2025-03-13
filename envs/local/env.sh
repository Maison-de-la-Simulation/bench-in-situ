SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
WORKING_DIR=${SCRIPT_DIR}/working_dir_$(date +%s)
SCHEFILE=${WORKING_DIR}/scheduler.json
BUILD_DIR=${WORKING_DIR}/sim/build
SIMULATION_BIN=${BUILD_DIR}/main
PYTHON_ENV=deisa
PYTHON_ENV_TARGZ=${SCRIPT_DIR}/../../output/
PYTHON_WORKING_DIR="output/tmp"
PYTHON_VERSION=${PYTHON_VERSION:-"310"}
PDI_BUILD_DIR=${WORKING_DIR}/pdi/build
PDI_INSTALL_DIR=${WORKING_DIR}/pdi/install
DASK_DISTRIBUTED__COMM__UCX__INFINIBAND="False"
DASK_NB_WORKERS=2
DASK_NB_THREAD_PER_WORKER=1 # if set to 0, will auto select depending on host (ncores)
DASK_WORKER_LOCAL_DIRECTORY="/tmp"
DASK_WORKER_PORT=8789

RED="\e[31m"
GREEN="\e[32m"
ENDCOLOR="\e[0m"

mkdir -p ${WORKING_DIR}
mkdir -p ${WORKING_DIR}/${PYTHON_ENV}
mkdir -p ${BUILD_DIR}
mkdir -p ${PDI_BUILD_DIR}

print_env() {
  echo "============="
  echo "SCRIPT_DIR=${SCRIPT_DIR}"
  echo "WORKING_DIR=${WORKING_DIR}"
  echo "SCHEFILE=${SCHEFILE}"
  echo "BUILD_DIR=${BUILD_DIR}"
  echo "PYTHON_ENV=${PYTHON_ENV}"
  echo "PYTHON_ENV_TARGZ=${PYTHON_ENV_TARGZ}"
  echo "PYTHON_VERSION=${PYTHON_VERSION}"
  echo "DASK_DISTRIBUTED__COMM__UCX__INFINIBAND=${DASK_DISTRIBUTED__COMM__UCX__INFINIBAND}"
  echo "DASK_NB_WORKERS=${DASK_NB_WORKERS}"
  echo "DASK_NB_THREAD_PER_WORKER=${DASK_NB_THREAD_PER_WORKER}"
  echo "============="
}

#expand aliases defined in the shell ~/.profile and ~/.bashrc
shopt -s expand_aliases
[ -f ~/.profile ] && source ~/.profile
[ -f ~/.bashrc ] && source ~/.bashrc
