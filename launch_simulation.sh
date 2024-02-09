#!/bin/bash

source env.sh
print_env

cd ${WORKING_DIR}

export PYTHONPATH=${WORKING_DIR}/${PYTHON_ENV}/lib/python3.11/site-packages
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/.local/lib/

echo "PYTHONPATH=${PYTHONPATH}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"

pdirun mpirun -np 1 ${BUILD_DIR}/main ../setup.ini ../io_deisa.yml

cd --
