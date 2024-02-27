SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
WORKING_DIR=${SCRIPT_DIR}/../working_dir
SCHEFILE=${WORKING_DIR}/scheduler.json
BUILD_DIR=${WORKING_DIR}/build
SIMULATION_BIN=${BUILD_DIR}/main
PYTHON_ENV=deisa
PYTHON_DEPS="deps"
PYTHON_VERSION=${PYTHON_VERSION:-"311"}
PDI_BUILD_DIR=${WORKING_DIR}/build/pdi
PDI_INSTALL_DIR=${WORKING_DIR}/pdi

mkdir -p ${WORKING_DIR}
mkdir -p ${WORKING_DIR}/${PYTHON_ENV}
mkdir -p ${PDI_BUILD_DIR}

print_env() {
	echo "============="
	echo "SCRIPT_DIR=${SCRIPT_DIR}"
	echo "WORKING_DIR=${WORKING_DIR}"
	echo "SCHEFILE=${SCHEFILE}"
	echo "BUILD_DIR=${BUILD_DIR}"
	echo "PYTHON_ENV=${PYTHON_ENV}"
	echo "PYTHON_ENV=${PYTHON_VERSION}"
	echo "============="
}
