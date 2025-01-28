#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env

cd ${WORKING_DIR}


echo "#######################"
echo "# Old Deisa"
echo "#######################"

# setup python environment
python3 -m venv ${PYTHON_ENV}

# activate python environment
source ${PYTHON_ENV}/bin/activate
pip install --upgrade pip

# untar and install dependencies
mkdir tmp_untar
tar xvf deisa_deps_py${PYTHON_VERSION}.tar.gz -C tmp_untar
pip install tmp_untar/*.whl
rm -rf tmp_untar

# install Deisa
pip install --no-index --no-build-isolation --no-deps ../lib/deisa

echo "TESTING DEPENDENCIES"
echo "Deisa.."
python3 -c "import deisa"
echo "HDF5.."
python3 -c "import h5py"
echo "Numpy.."
python3 -c "import numpy"
echo "Dask.."
python3 -c "import dask"
echo "Distributed.."
python3 -c "import distributed"

deactivate


echo "#######################"
echo "# New Deisa"
echo "#######################"
PYTHON_ENV="${PYTHON_ENV}_new"

# setup python environment
python3 -m venv ${PYTHON_ENV}

# activate python environment
source ${PYTHON_ENV}/bin/activate
pip install --upgrade pip

# untar and install dependencies
mkdir tmp_untar
tar xvf new_deisa_deps_py${PYTHON_VERSION}.tar.gz -C tmp_untar
pip install tmp_untar/*.whl
rm -rf tmp_untar

# install Deisa
pip install --no-index --no-build-isolation --no-deps ../lib/new_deisa

echo "TESTING DEPENDENCIES"
echo "Deisa.."
python3 -c "import deisa"
echo "HDF5.."
python3 -c "import h5py"
echo "Numpy.."
python3 -c "import numpy"
echo "Dask.."
python3 -c "import dask"
echo "Distributed.."
python3 -c "import distributed"

deactivate

cd --
