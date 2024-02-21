#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
cd ${SCRIPT_DIR}
source env.sh

./offline_python_env_setup.sh
./build_pdi.sh
./build_simulation.sh


