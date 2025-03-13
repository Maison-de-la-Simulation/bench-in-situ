#!/bin/bash

############################
# This script:
# - builds python virtual environment (new and old deisa)
# - build and install PDI
# - build the simulation
############################

set -xeu
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

source ${SCRIPT_DIR}/env.sh
cd ${WORKING_DIR}
print_env

python3 --version
type python3

echo "###########################"
echo "# Building OLD Deisa venv #"
echo "###########################"

echo "Setting up python venv"
python3 -m venv ${PYTHON_ENV}
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "Python venv failed !"
  exit 1
fi

sync

echo "Activating python venv"
source ${PYTHON_ENV}/bin/activate
if [ -z "${VIRTUAL_ENV}" ]; then
  echo "Could not activate python environment !"
  exit 1
fi


pip install --upgrade pip
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "Python pip upgrade failed !"
  exit 1
fi

DEPS_FILE="${PYTHON_ENV_TARGZ}/deisa_deps_py${PYTHON_VERSION}.tar.gz"
if [ ! -f "$DEPS_FILE" ]; then
    echo "File not found: $DEPS_FILE"
    exit 1
fi

# untar and install dependencies
mkdir tmp_untar
tar xvf ${DEPS_FILE} -C tmp_untar
pip install tmp_untar/*.whl
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "Python dependency install failed !"
  exit 1
fi
rm -rf tmp_untar

# install Deisa
pip install --no-index --no-build-isolation --no-deps ${SCRIPT_DIR}/../../lib/deisa
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "Python deisa install failed !"
  exit 1
fi

deactivate
echo "venv build successfully"

echo "###########################"
echo "# Building NEW Deisa venv #"
echo "###########################"

PYTHON_ENV="${PYTHON_ENV}_new"

# setup python environment
python3 -m venv ${PYTHON_ENV}
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "Python venv failed for NEW Deisa !"
  exit 1
fi

# activate python environment
source ${PYTHON_ENV}/bin/activate
if [ -z "${VIRTUAL_ENV}" ]; then
  echo "Could not activate python environment !"
  exit 1
fi

pip install --upgrade pip
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "Python pip upgrade failed !"
  exit 1
fi

DEPS_FILE="${PYTHON_ENV_TARGZ}/new_deisa_deps_py${PYTHON_VERSION}.tar.gz"
if [ ! -f "$DEPS_FILE" ]; then
    echo "File not found: $DEPS_FILE"
    exit 1
fi

# untar and install dependencies
mkdir tmp_untar
tar xvf ${DEPS_FILE} -C tmp_untar
pip install tmp_untar/*.whl
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "Python dependency install failed !"
  exit 1
fi
rm -rf tmp_untar

# install Deisa
pip install --no-index --no-build-isolation --no-deps ${SCRIPT_DIR}/../../lib/new_deisa
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "Python deisa install failed !"
  exit 1
fi

deactivate
echo "venv build successfully"

echo "################"
echo "# Building PDI #"
echo "################"

echo "PYTHONPATH=${PYTHONPATH}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
echo "LD_PRELOAD=${LD_PRELOAD}"
echo "PYTHON_ENV=${PYTHON_ENV}"

# TODO: what python env should be used ? new or old ?
source ${PYTHON_ENV}/bin/activate

export LD_LIBRARY_PATH=${PYTHON_ENV}/lib:${LD_LIBRARY_PATH}

python --version

cd ${PDI_BUILD_DIR}

cmake -DCMAKE_INSTALL_PREFIX=${PDI_INSTALL_DIR} \
  -DBUILD_DECL_HDF5_PLUGIN=OFF -DBUILD_DECL_NETCDF_PLUGIN=OFF \
  -DBUILD_HDF5_PARALLEL=OFF -DBUILD_SHARED_LIBS=ON -DBUILD_FORTRAN=OFF \
  -DBUILD_NETCDF_PARALLEL=OFF \
  -DUSE_yaml=EMBEDDED -DUSE_pybind11=EMBEDDED -DUSE_paraconf=EMBEDDED \
  -DBUILD_BENCHMARKING=OFF -DBUILD_TESTING=OFF \
  -DBUILD_SET_VALUE_PLUGIN=OFF -DBUILD_DECL_NETCDF_PLUGIN=OFF -DBUILD_USER_CODE_PLUGIN=OFF \
  -DBUILD_PYTHON=ON -DBUILD_DEISA_PLUGIN=ON \
  ${SCRIPT_DIR}/../../lib/pdi
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "PDI cmake failed !"
  exit 1
fi

make -j $(nproc)
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "PDI make failed !"
  exit 1
fi
make install

deactivate
cd ${WORKING_DIR}

echo "PDI installed"

echo "#######################"
echo "# Building SIMULATION #"
echo "#######################"

source ${PDI_INSTALL_DIR}/share/pdi/env.sh

cd ${BUILD_DIR}


#cmake -DCMAKE_BUILD_TYPE=Release -DSESSION=MPI_SESSION -DKokkos_ENABLE_OPENMP=ON -DEuler_ENABLE_PDI=ON ../..
#    -DCMAKE_CXX_STANDARD=17 \

cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DEuler_ENABLE_PDI=ON \
    -DKokkos_ENABLE_OPENMP=ON \
    -DKokkos_ENABLE_SERIAL=OFF \
    -DKokkos_ENABLE_CUDA=OFF \
    -DSESSION=MPI_SESSION \
    ${SCRIPT_DIR}/../../simulation
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "Simulation cmake failed !"
  exit 1
fi

make -j $(nproc) 
rc=$?
if [ ${rc} -ne 0 ]; then
  echo "Simulation make failed !"
  exit 1
fi

cd ${WORKING_DIR}

echo "Simulation built successfully"
echo "DONE !"
