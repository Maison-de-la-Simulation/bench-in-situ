#!/bin/bash

source env.sh
print_env

exit

#SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
#WORKING_DIR=${SCRIPT_DIR}/working_dir # TODO: Add datetime ?
#SCHEFILE=${WORKING_DIR}/scheduler.json
#BUILD_DIR=${SCRIPT_DIR}/out/Debug

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

# Simulation uses Deisa, which requires a properly configured Python environment.
PYTHONPATH=./python_env
mkdir -p ${PYTHONPATH}

PYTHONPATH=${PYTHONPATH} ../${BUILD_DIR}/main

cd --
