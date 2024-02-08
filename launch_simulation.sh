#!/bin/bash

source env.sh
print_env

cd ${WORKING_DIR}

PYTHONPATH=${PYTHON_ENV} ${BUILD_DIR}/main setup.ini io_deisa.yml

cd --
