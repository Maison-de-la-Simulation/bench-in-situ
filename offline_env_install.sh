#!/bin/bash

source env.sh
print_env

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

# setup python environment
python -m venv ${PYTHON_ENV}

# activate python environment
source ${PYTHON_ENV}/bin/activate

pip install deps/*.whl
pip install --no-index --no-build-isolation deps/*.zip

# install Deisa
pip install --no-index --no-build-isolation --no-deps ../lib/deisa

cd --