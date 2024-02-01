mkdir benchInSitu
~/miniconda3/bin/conda create python=3.7 --prefix ./condaEnv
conda activate ./condaEnv

git clone --recurse-submodules https://github.com/Maison-de-la-Simulation/bench-in-situ.git
cd bench-in-situ/lib/pdi

git checkout deisa_plugin
sed -i 's/project(pdi_deisa_plugin LANGUAGES CXX)/project(pdi_deisa_plugin LANGUAGES C CXX)/' plugins/deisa/CMakeLists.txt

mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../../install_pdi -DUSE_HDF5=EMBEDDED -DBUILD_HDF5_PARALLEL=ON  -DUSE_yaml=EMBEDDED -DUSE_paraconf=EMBEDDED -DBUILD_SHARED_LIBS=ON -DBUILD_FORTRAN=OFF -DBUILD_BENCHMARKING=OFF -DBUILD_SET_VALUE_PLUGIN=OFF -DBUILD_TESTING=OFF -DBUILD_DECL_NETCDF_PLUGIN=OFF -DBUILD_USER_CODE_PLUGIN=ON -DBUILD_DEISA_PLUGIN=ON -DBUILD_PYTHON=ON ..
make -j 8
make install

cd ../../../
mkdir build && cd build
cmake -DSESSION=MPI_SESSION -DKokkos_ENABLE_OPENMP=ON -DEuler_ENABLE_PDI=ON ..
make -j 4