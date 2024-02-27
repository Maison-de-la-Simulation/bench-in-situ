#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
cd ${SCRIPT_DIR}
source env.sh

./offline_python_env_setup.sh
./build_pdi.sh
./build_simulation.sh

# copy files to working directory
cd ..
cp deisa_deps_py* ${WORKING_DIR}
cp -r in-situ ${WORKING_DIR}
#cp io.yml io_chkpt.yml io_deisa.yml ${WORKING_DIR}

