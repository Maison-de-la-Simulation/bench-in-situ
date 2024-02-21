#!/bin/bash

source env.sh
print_env
source ${WORKING_DIR}/install_pdi/share/pdi/env.sh

cd ${WORKING_DIR}
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

cmake -CMAKE_BUILD_TYPE=Release -DSESSION=MPI_SESSION -DKokkos_ENABLE_OPENMP=ON -DEuler_ENABLE_PDI=ON ../..

make -j8

cd --

