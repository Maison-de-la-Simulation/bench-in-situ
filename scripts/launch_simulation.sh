#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env
source ${WORKING_DIR}/${PYTHON_ENV}/bin/activate
source ${PDI_INSTALL_DIR}/share/pdi/env.sh


export PYTHONPATH=${WORKING_DIR}/${PYTHON_ENV}/lib/python3.11/site-packages:${PYTHONPATH}
#export LD_LIBRARY_PATH=${PDI_INSTALL_DIR}/lib:${WORKING_DIR}/${PYTHON_ENV}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${WORKING_DIR}/${PYTHON_ENV}/lib:${LD_LIBRARY_PATH}
#export LD_PRELOAD=${PDI_INSTALL_DIR}/lib64/libyaml.so # note: only for mac

echo "PYTHONPATH=${PYTHONPATH}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
#echo "LD_PRELOAD=${LD_PRELOAD}"

cp ../io_chkpt.yml ../io_deisa.yml ${BUILD_DIR}
cp ../setup.ini ${BUILD_DIR}
cp ${SCHEFILE} ${BUILD_DIR}

cd ${BUILD_DIR}

./main setup.ini io_deisa.yml

#pdirun mpirun -np 1 ${BUILD_DIR}/main ../setup.ini ../io_deisa.yml
cd --
