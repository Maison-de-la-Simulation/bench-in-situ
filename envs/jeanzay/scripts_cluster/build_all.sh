#!/bin/bash
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
cd ${SCRIPT_DIR}
source env.sh

if [ "${CLUSTER_NAME}" == "jeanzay" ]; then
    source ../modules_${GPU_ARCH}.env
fi

./build_pdi.sh
./build_simulation.sh

# copy files to working directory
cd ..
cp -r ${BENCH_ROOT_DIR}/in-situ ${WORKING_DIR}
#cp io.yml io_chkpt.yml io_deisa.yml ${WORKING_DIR}

