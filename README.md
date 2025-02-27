# bench-in-situ

This benchmark measures the impact of insitu analytics on the I/O capability of a system.
While the HPC simulation is running, periodic HDF5 checkpoints are performed, generating large files for each simulation time-step.
At the same time, insitu analytics is computed and multiple small files are produced.





## To clone the project

``` bash
git clone --recurse-submodules https://github.com/Maison-de-la-Simulation/bench-in-situ.git
```

or

``` bash
git clone https://github.com/Maison-de-la-Simulation/bench-in-situ.git 
git submodule init
git submodule update
```

### Modules needed (on Ruche)

``` bash
module load gcc/11.2.0/gcc-4.8.5
module load openmpi/4.1.1/gcc-11.2.0
module load cuda/11.7.0/gcc-11.2.0
module load cmake/3.21.4/gcc-11.2.0 
```

### Modules needed (on Jean-Zay)

``` bash
module load cmake/3.21.3
module load intel-compilers/19.0.4
module load intel-mpi/2019.4
```


### To build PDI

``` bash
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../../install_pdi -DUSE_HDF5=EMBEDDED -DBUILD_HDF5_PARALLEL=ON  -DUSE_yaml=EMBEDDED -DUSE_paraconf=EMBEDDED -DBUILD_SHARED_LIBS=ON -DBUILD_FORTRAN=OFF -DBUILD_BENCHMARKING=OFF -DBUILD_SET_VALUE_PLUGIN=OFF -DBUILD_TESTING=OFF -DBUILD_DECL_NETCDF_PLUGIN=OFF -DBUILD_USER_CODE_PLUGIN=ON ..
make -j 8
make install
source $PWD/../../install_pdi/share/pdi/env.sh
```

### To build bench

``` bash
mkdir build && cd build
cmake -DSESSION=MPI_SESSION -DKokkos_ENABLE_OPENMP=ON -DEuler_ENABLE_PDI=ON ..
make -j 4
```

### To run bench

``` bash
./main ../setup.ini ../io_chkpt.yml
```

### Offline Installation

* Download the Python environment on the online machine : `./scripts/offline_python_env_download.sh`
* Copy the whole repo on the offline machine
* Source the env file and ruche env file : `source scripts/env.sh && source envs/PLATFORM/modules.env`
* unzip the `deisa_.tar.gz` file according to the Python version used on the offline machine : `tar -xzvf working_dir/deisa_deps_py${PYTHON_VERSION}.tar.gz -C ${WORKING_DIR}`
* Install the Python environment on the offline machine : `./scripts/offline_python_env_setup.sh`
* Install PDI : `./scripts/build_pdi.sh `
* Install ARK-MHD : `./scripts/build_simulation.sh `


# Using Guix

``` bash
guix shell --pure python coreutils grep glibc zlib tar gzip bash lesspipe sed guix
export LD_PRELOAD=$GUIX_ENVIRONMENT/lib/libz.so
```



