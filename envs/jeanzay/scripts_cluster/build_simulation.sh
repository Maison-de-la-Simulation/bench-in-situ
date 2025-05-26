#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env
echo  "PDI_INSTALL_DIR=${PDI_INSTALL_DIR}"
source ${PDI_INSTALL_DIR}/share/pdi/env.sh

# If we want only to recompile the simulation code.
if [ "${CLUSTER_NAME}" == "jeanzay" ]; then
    if [ "${GPU_ARCH}" == "V100" ] ||  [ "${GPU_ARCH}" == "A100" ] ||  [ "${GPU_ARCH}" == "H100" ]; then
        source ../modules_${GPU_ARCH}.env
    fi
fi


cd ${WORKING_DIR}
#source ${PYTHON_ENV}/bin/activate

mkdir -p ${BUILD_SIMULATION_DIR}
cd ${BUILD_SIMULATION_DIR}

if [ "${CLUSTER_NAME}" == "jeanzay" ]; then

    if [ "${GPU_ARCH}" == "V100" ]; then
        # install v100
        # skx, avx512, avx2 need to be tested
        cmake \
            -DCMAKE_BUILD_TYPE=Release \
            -DEuler_ENABLE_PDI=ON \
            -DKokkos_ENABLE_OPENMP=ON \
            -DKokkos_ENABLE_SERIAL=OFF \
            -DKokkos_ENABLE_CUDA=ON \
            -DKokkos_ARCH_VOLTA70=ON \
            -DKokkos_ARCH_SKX=ON \
            -DKokkos_ARCH_AVX512XEON=OFF \
            -DKokkos_ARCH_AVX2=OFF \
            -DSESSION=MPI_SESSION \
            ../../../../../simulation
    elif [ "${GPU_ARCH}" == "A100" ]; then
        # install a100
        # zen3 and AMD_AVX2 need to be tested
        cmake \
            -DCMAKE_BUILD_TYPE=Release \
            -DEuler_ENABLE_PDI=ON \
            -DKokkos_ENABLE_OPENMP=ON \
            -DKokkos_ENABLE_SERIAL=OFF \
            -DKokkos_ENABLE_CUDA=ON \
            -DKokkos_ARCH_AMPERE80=ON \
            -DKokkos_ARCH_ZEN3=ON \
            -DKokkos_ARCH_AMD_AVX2=ON \
            -DSESSION=MPI_SESSION \
            ../../../../../simulation

    elif [ "${GPU_ARCH}" == "H100" ]; then
        # install h100
        # SPR,AVX2, AVX512XEON need to be tested
        cmake \
            -DCMAKE_BUILD_TYPE=Release \
            -DEuler_ENABLE_PDI=ON \
            -DKokkos_ENABLE_OPENMP=ON \
            -DKokkos_ENABLE_SERIAL=OFF \
            -DKokkos_ENABLE_CUDA=ON \
            -DKokkos_ARCH_HOPPER90=ON \
            -DKokkos_ARCH_SPR=ON \
            -DKokkos_ARCH_AVX512XEON=ON \
            -DKokkos_ARCH_AVX2=ON \
            -DSESSION=MPI_SESSION \
            ../../../../../simulation
    else
        echo "jeanzay is only with gpu"
        exit 1
    fi

else
    # install cpu
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DEuler_ENABLE_PDI=ON \
        -DKokkos_ENABLE_OPENMP=ON \
        -DKokkos_ENABLE_SERIAL=OFF \
        -DKokkos_ENABLE_CUDA=OFF \
        -DKokkos_ARCH_AMPERE80=OFF \
        -DKokkos_ARCH_PASCAL60=OFF \
        -DKokkos_ARCH_ZEN3=OFF \
        -DKokkos_ENABLE_HIP=OFF \
        -DKokkos_ARCH_VEGA90A=OFF \
        -DSESSION=MPI_SESSION \
        ../..
fi

#Ampere : ruche a100
#Pascal : ruche p100
#Zen3, HIP and Vega : Adastra
make -j 8 #$(nproc)

#deactivate

cd --

