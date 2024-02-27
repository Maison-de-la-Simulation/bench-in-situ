#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env

cd ${WORKING_DIR}

# setup python environment
python -m venv ${PYTHON_ENV}

# activate python environment
source ${PYTHON_ENV}/bin/activate


tar xvf deisa_deps_py${PYTHON_VERSION}.tar.gz
pip install ${PYTHON_DEPS}/*.whl

# install Deisa
pip install --no-index --no-build-isolation --no-deps ../lib/deisa

echo "TESTING DEPENDENCIES"
echo "Deisa.."
python -c "import deisa"
echo "Numpy.."
python -c "import numpy"
echo "Dask.."
python -c "import dask"
echo "Distributed.."
python -c "import distributed"

cd --
