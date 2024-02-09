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

# install Deisa
pip install --no-index --no-build-isolation --no-deps ../lib/deisa

echo "TESTING DEPENDENCIES"
echo "Deisa :"
python -m import deisa
echo "Numpy :"
python -m import numpy
echo "Dask :"
python -m import Dask
echo "Distributed :"
python -m import distributed
cd --