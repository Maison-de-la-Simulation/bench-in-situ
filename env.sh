SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
WORKING_DIR=${SCRIPT_DIR}/working_dir # TODO: Add datetime ?
SCHEFILE=${WORKING_DIR}/scheduler.json
BUILD_DIR=${SCRIPT_DIR}/out/Debug
PYTHON_ENV=deisa
SIMULATION_BIN=${BUILD_DIR}/main
PYTHON_DEPS=${WORKING_DIR}/deps

mkdir -p ${WORKING_DIR}
mkdir -p ${WORKING_DIR}/${PYTHON_ENV}

print_env() {
	echo "============="
	echo "SCRIPT_DIR=${SCRIPT_DIR}"
	echo "WORKING_DIR=${WORKING_DIR}"
	echo "SCHEFILE=${SCHEFILE}"
	echo "BUILD_DIR=${BUILD_DIR}"
	echo "PYTHON_ENV=${PYTHON_ENV}"
	echo "============="
}
