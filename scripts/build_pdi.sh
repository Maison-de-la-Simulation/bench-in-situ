#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env

cd ${WORKING_DIR}
source ${PYTHON_ENV}/bin/activate

python --version

cd ${PDI_BUILD_DIR}

cmake -DCMAKE_INSTALL_PREFIX=${PDI_INSTALL_DIR} \
  -DUSE_HDF5=EMBEDDED -DUSE_yaml=EMBEDDED -DUSE_pybind11=EMBEDDED -DUSE_paraconf=EMBEDDED \
  -DBUILD_HDF5_PARALLEL=OFF -DBUILD_SHARED_LIBS=ON -DBUILD_FORTRAN=OFF \
  -DBUILD_BENCHMARKING=OFF -DBUILD_TESTING=OFF \
  -DBUILD_SET_VALUE_PLUGIN=OFF -DBUILD_DECL_NETCDF_PLUGIN=OFF -DBUILD_USER_CODE_PLUGIN=ON -DBUILD_PYTHON=ON -DBUILD_DEISA_PLUGIN=ON \
  ../../../lib/pdi

make -j8
make install

source ${PDI_INSTALL_DIR}/share/pdi/env.sh

cd --

