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

## Modules needed
### Modules needed on Ruche
``` bash
module load gcc/11.2.0/gcc-4.8.5
module load openmpi/4.1.1/gcc-11.2.0
module load cuda/11.7.0/gcc-11.2.0
module load cmake/3.21.4/gcc-11.2.0 
```

### Modules needed on Jean-Zay
``` bash
module load cmake/3.21.3
module load intel-compilers/19.0.4
module load intel-mpi/2019.4
```

### Modules needed on Adastra for MI250X
```
module load PrgEnv-cray
module load amd-mixed/6.1.2
module load rocm/6.1.2
module load craype-accel-amd-gfx90a
module load craype-x86-trento
module load cray-python
```

### Modules needed on Adastra for MI300A
```
module load PrgEnv-cray
module load amd-mixed/6.1.2
module load rocm/6.1.2
module load craype-accel-amd-gfx942
module load craype-x86-genoa
module load cray-python
```

## To build PDI
From `bench-in-situ/lib/pdi`
``` bash
mkdir build && cd build
cmake \
	-DCMAKE_INSTALL_PREFIX=$PWD/../../install_pdi \
	-DUSE_HDF5=EMBEDDED \
	-DBUILD_HDF5_PARALLEL=ON \
	-DUSE_yaml=EMBEDDED \
	-DUSE_paraconf=EMBEDDED \
	-DBUILD_SHARED_LIBS=ON \
	-DBUILD_FORTRAN=OFF \
	-DBUILD_BENCHMARKING=OFF \
	-DBUILD_SET_VALUE_PLUGIN=OFF \
	-DBUILD_TESTING=OFF \
	-DBUILD_DECL_NETCDF_PLUGIN=OFF \
	-DBUILD_USER_CODE_PLUGIN=ON \
	..
make -j8
make install
source $PWD/../../install_pdi/share/pdi/env.sh
```

## To build bench
From `bench-in-situ/simulation`
``` bash
mkdir build && cd build
cmake \
	-DSESSION=MPI_SESSION \
	-DKokkos_ENABLE_OPENMP=ON \
	-DEuler_ENABLE_PDI=ON ..
make -j8
```

### To build bench on Adastra
For MI250, use the following from `bench-in-situ/simulation`:
``` bash
export ${ACTIVE_PROJECT}
source $PWD/../lib/install_pdi/share/pdi/env.sh
source $PWD/../envs/adastra/modules.env
source $PWD/../envs/adastra/modulesMI250.env
mkdir build && cd build
cmake \
	-DSESSION=MPI_SESSION \
	-DKokkos_ARCH_AMD_GFX90A=ON \
	-DKokkos_ENABLE_HIP=ON \
	-DKokkos_ENABLE_OPENMP=ON \
	-DEuler_ENABLE_PDI=ON \
	-DCMAKE_CXX_FLAGS="-fopenmp" \
	..
make -j16
```
For MI300, use the following from `bench-in-situ/simulation`:
``` bash
export ${ACTIVE_PROJECT}
source $PWD/../lib/install_pdi/share/pdi/env.sh
source $PWD/../envs/adastra/modules.env
source $PWD/../envs/adastra/modulesMI300.env
mkdir build && cd build
cmake \
	-DSESSION=MPI_SESSION \
	-DKokkos_ARCH_AMD_GFX942=ON \
	-DKokkos_ENABLE_HIP=ON \
	-DKokkos_ENABLE_OPENMP=ON \
	-DEuler_ENABLE_PDI=ON \
	-DCMAKE_CXX_FLAGS="-fopenmp" \
	..
make -j16
```

### To build bench on Jean-Zay
Use the following from `bench-in-situ/simulation`:
``` bash
mkdir build && cd build
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DEuler_ENABLE_MPI_CUDA_AWARE=ON \
    -DEuler_ENABLE_PDI=ON \
    -DKokkos_ENABLE_OPENMP=ON \
    -DKokkos_ENABLE_SERIAL=OFF \
    -DKokkos_ENABLE_CUDA=ON \
    -DKokkos_ARCH_AMPERE80=OFF \
    -DKokkos_ARCH_PASCAL60=OFF \
    -DKokkos_ARCH_VOLTA70=ON \
    -DKokkos_ARCH_ZEN3=OFF \
    -DKokkos_ENABLE_HIP=OFF \
    -DKokkos_ARCH_VEGA90A=OFF \
    -DSESSION=MPI_SESSION \
    -DCMAKE_CXX_STANDARD=17 \
    ..
make -j8
```

## To run bench
Run manually from the folder `simulation` using:
``` bash
./main ../setup.ini ../io_chkpt.yml
```
Or use a script from a subdirectory of `envs` with sbatch (or bash), 
for example:
``` bash
bash ./envs/adastra/strong/launcherBench.sh MI250
```

## Offline Installation

* Download the Python environment on the online machine : `./scripts/offline_python_env_download.sh`
* Copy the whole repository on the offline machine
* Source the `.env` file and ruche `.env` file : `source scripts/env.sh && source envs/PLATFORM/modules.env`
* Unzip the `deisa_.tar.gz` file according to the Python version used on the offline machine : `tar -xzvf working_dir/deisa_deps_py${PYTHON_VERSION}.tar.gz -C ${WORKING_DIR}`
* Install the Python environment on the offline machine : `./scripts/offline_python_env_setup.sh`
* Install PDI : `./scripts/build_pdi.sh `
* Install ARK-MHD : `./scripts/build_simulation.sh `


# Using Guix

``` bash
guix shell --pure python coreutils grep glibc zlib tar gzip bash lesspipe sed guix
export LD_PRELOAD=$GUIX_ENVIRONMENT/lib/libz.so
```



