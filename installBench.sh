# mkdir benchInSitu
# ~/miniconda3/bin/conda create python=3.7 --prefix ./condaEnv
# conda activate ./condaEnv

cd bench-in-situ/lib/pdi

#git checkout deisa_plugin
#sed -i 's/project(pdi_deisa_plugin LANGUAGES CXX)/project(pdi_deisa_plugin LANGUAGES C CXX)/' plugins/deisa/CMakeLists.txt

mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../../install_pdi -DUSE_DEFAULT=EMBEDDED -DBUILD_FORTRAN=OFF -DBUILD_BENCHMARKING=OFF -DBUILD_SET_VALUE_PLUGIN=OFF -DBUILD_TESTING=OFF -DBUILD_USER_CODE_PLUGIN=ON -DBUILD_DEISA_PLUGIN=ON -DBUILD_PYTHON=ON ..
make -j 8
make install
source ../../install_pdi/share/pdi/env.sh

cd ../../../
mkdir build && cd build
cmake -DSESSION=MPI_SESSION -DKokkos_ENABLE_OPENMP=ON -DEuler_ENABLE_PDI=ON ..
make -j 4

