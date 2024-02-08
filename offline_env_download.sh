#!/bin/bash

source env.sh
print_env

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

# TODO: create conda python 3.9
# 1) create conda env with python 3.9: conda create -y --name py39 python=3.9.10
# 2) activate environment
# 3) pip download
# 4) tar gz -> deisa_deps_py39.tar.gz
# 5) cleanup: deactivate conda env, rm conda env, downloaded files

# TODO: create conda python 3.10
# TODO: create conda python 3.11

# setup python environment
python -m venv ${PYTHON_ENV}

# activate python environment
source ${PYTHON_ENV}/bin/activate

# install Deisa
pip install ../lib/deisa

# TODO: offline setup

pip download ../lib/deisa/ -d deps/
pip download wheel versioneer -d deps/

deactivate
rm -r ${PYTHON_ENV}
cd --

