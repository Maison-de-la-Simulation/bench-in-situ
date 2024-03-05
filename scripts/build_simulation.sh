#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env
source ${PDI_INSTALL_DIR}/share/pdi/env.sh

cd ${WORKING_DIR}
source ${PYTHON_ENV}/bin/activate

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}


#cmake -DCMAKE_BUILD_TYPE=Release -DSESSION=MPI_SESSION -DKokkos_ENABLE_OPENMP=ON -DEuler_ENABLE_PDI=ON ../..

cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=17 \
    -DEuler_ENABLE_PDI=ON \
    -DKokkos_ENABLE_SERIAL=ON \
    -DKokkos_ENABLE_CUDA=OFF \
    -DKokkos_ARCH_AMPERE80=OFF \
    -DKokkos_ARCH_PASCAL60=OFF \
    -DKokkos_ARCH_ZEN3=ON \
    -DKokkos_ENABLE_HIP=ON \
    -DKokkos_ARCH_VEGA90A=ON \
    -DSESSION=MPI_SESSION \
    ../..
#Ampere : ruche a100
#Pascal : ruche p100
#Zen3, HIP and Vega : Adastra
make -j $(nproc) 

cd --

