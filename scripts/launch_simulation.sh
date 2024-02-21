#!/bin/bash

source env.sh
print_env
source lib/install_pdi/share/pdi/env.sh


export PYTHONPATH=${WORKING_DIR}/${PYTHON_ENV}/lib/python3.11/site-packages
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/.local/lib/

echo "PYTHONPATH=${PYTHONPATH}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"

cp io_chkpt.yml ${BUILD_DIR}
cp setup.ini ${BUILD_DIR}
cd ${BUILD_DIR}

export LD_PRELOAD=${WORKING_DIR}/install_pdi/lib64/libyaml.so

./main ./setup.ini ./io_chkpt.yml

#pdirun mpirun -np 1 ${BUILD_DIR}/main ../setup.ini ../io_deisa.yml
cd --
