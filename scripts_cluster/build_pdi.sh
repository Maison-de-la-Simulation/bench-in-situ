#!/bin/bash
for i in "$@"
do
case $i in
    -gpu=*)
    GPU_ARCH="${i#*=}"
    ;;
    -cluster=*)
    CLUSTER_NAME="${i#*=}"
    ;;
    -deisa=*)
    PYTHON_DEISA="${i#*=}"
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
    # unknown option
    ;;
esac
done

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh -gpu=${GPU_ARCH} -cluster=${CLUSTER_NAME}
print_env

cd ${WORKING_DIR}

##source ${PYTHON_ENV}/bin/activate

if [ "${PYTHON_DEISA}" == "ON" ]; then
  source ${PYTHON_ENV}/bin/activate
  export LD_LIBRARY_PATH=${PYTHON_ENV}/lib:${LD_LIBRARY_PATH}
fi

python3 --version

cd ${PDI_BUILD_DIR}

cmake -DCMAKE_INSTALL_PREFIX=${PDI_INSTALL_DIR} \
  -DUSE_HDF5=EMBEDDED -DUSE_yaml=EMBEDDED -DUSE_pybind11=EMBEDDED -DUSE_paraconf=EMBEDDED \
  -DBUILD_HDF5_PARALLEL=ON -DBUILD_SHARED_LIBS=ON -DBUILD_FORTRAN=OFF \
  -DBUILD_BENCHMARKING=OFF -DBUILD_TESTING=OFF \
  -DBUILD_SET_VALUE_PLUGIN=ON -DBUILD_DECL_NETCDF_PLUGIN=OFF -DBUILD_USER_CODE_PLUGIN=ON -DBUILD_PYTHON=ON -DBUILD_DEISA_PLUGIN=ON \
  ../../../lib/pdi

make -j 8 #$(nproc)
make install

source ${PDI_INSTALL_DIR}/share/pdi/env.sh

if [ "${PYTHON_DEISA}" == "ON" ]; then
  deactivate
fi

cd --

