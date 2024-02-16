#!/bin/bash

source env.sh
print_env
source ${WORKING_DIR}/install_pdi/share/pdi/env.sh

cd ${WORKING_DIR}
mkdir -p build
cd build

cmake -DSESSION=MPI_SESSION -DKokkos_ENABLE_OPENMP=ON -DEuler_ENABLE_PDI=ON ../..

make -j8


cd --

