# bench-in-situ
## To clone the project

```
git clone --recurse-submodules https://github.com/Maison-de-la-Simulation/bench-in-situ.git
```

or

```
git clone https://github.com/Maison-de-la-Simulation/bench-in-situ.git 
git submodule init
git submodule update
```
### To build PDI
```
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../../install_pdi -DUSE_HDF5=EMBEDDED -DBUILD_HDF5_PARALLEL=ON  -DUSE_yaml=EMBEDDED -DUSE_paraconf=EMBEDDED -DBUILD_SHARED_LIBS=ON -DBUILD_FORTRAN=OFF -DBUILD_BENCHMARKING=OFF -DBUILD_SET_VALUE_PLUGIN=OFF -DBUILD_TESTING=OFF -DBUILD_DECL_NETCDF_PLUGIN=OFF -DBUILD_USER_CODE_PLUGIN=ON ..
make -j 8
make install
```

### To build bench
```
mkdir build && cd build
cmake -DSESSION=MPI_SESSION -DKokkos_ENABLE_OPENMP=ON -DEuler_ENABLE_PDI=ON ..
make -j 4
```
