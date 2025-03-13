#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

SCHEFILE=scheduler.json
PREFIX=bench_insitu

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores

cd ${WORKING_DIR}

source ${SCRIPT_DIR}/env.sh
source ${PDI_INSTALL_DIR}/share/pdi/env.sh

run_dask () {
  dask scheduler --scheduler-file=$SCHEFILE &
  sleep 3
  sync
  
  echo "Starting ${DASK_NB_WORKERS} Dask workers with ${DASK_NB_THREAD_PER_WORKER} threads per worker."
  # --nworkers sets the number of worker processes
  # --nthreads sets the number of threads per worker process
  dask worker --nworkers ${DASK_NB_WORKERS} --nthreads ${DASK_NB_THREAD_PER_WORKER}\
    --local-directory ${DASK_WORKER_LOCAL_DIRECTORY}\
    --scheduler-file=$SCHEFILE &
}


if [ "$1" = "old" ]; then
  echo "#####################"
  echo "# Running OLD Deisa #"
  echo "#####################"
  echo "Setting up Dask environement"

  # setup python environment
  python3 -m venv ${PYTHON_ENV}
  
  # activate python environment
  source ${PYTHON_ENV}/bin/activate
  pip install --upgrade pip
  
  DEPS_FILE="${PYTHON_ENV_TARGZ}/deisa_deps_py${PYTHON_VERSION}.tar.gz"
  if [ ! -f "$DEPS_FILE" ]; then
      echo "File not found: $DEPS_FILE"
      exit 1
  fi
  
  # untar and install dependencies
  mkdir tmp_untar
  tar xvf ${DEPS_FILE} -C tmp_untar
  pip install tmp_untar/*.whl
  rm -rf tmp_untar
  
  # install Deisa
  pip install --no-index --no-build-isolation --no-deps ${SCRIPT_DIR}/../../lib/deisa

  run_dask
  python3 ${SCRIPT_DIR}/../../in-situ/fft.py $SCHEFILE &

elif [ "$1" = "new" ]; then
  echo "#####################"
  echo "# Running NEW Deisa #"
  echo "#####################"
  echo "Setting up Dask environement"

  PYTHON_ENV="${PYTHON_ENV}_new"

  # setup python environment
  python3 -m venv ${PYTHON_ENV}
  
  # activate python environment
  source ${PYTHON_ENV}/bin/activate
  pip install --upgrade pip
  
  DEPS_FILE="${PYTHON_ENV_TARGZ}/new_deisa_deps_py${PYTHON_VERSION}.tar.gz"
  if [ ! -f "$DEPS_FILE" ]; then
      echo "File not found: $DEPS_FILE"
      exit 1
  fi
  
  # untar and install dependencies
  mkdir tmp_untar
  tar xvf ${DEPS_FILE} -C tmp_untar
  pip install tmp_untar/*.whl
  rm -rf tmp_untar
  
  # install Deisa
  pip install --no-index --no-build-isolation --no-deps ${SCRIPT_DIR}/../../lib/new_deisa
    
  run_dask
  python3 ${SCRIPT_DIR}/../../in-situ/bench_deisa.py $SCHEFILE &
else
  echo "Unknown option. accepted values are: old, new"
  exit 1
fi



echo "PYTHONPATH=${PYTHONPATH}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
echo "LD_PRELOAD=${LD_PRELOAD}"

cp ${SCRIPT_DIR}/deisa/io.yml ${BUILD_DIR}
cp ${SCRIPT_DIR}/deisa/setup.ini ${BUILD_DIR}
cp ${SCHEFILE} ${BUILD_DIR}

cd ${BUILD_DIR}

echo "Running simulation"
#perf record mpirun -np 2 ./main setup.ini io_deisa.yml
mpirun -np 4 ${SIMULATION_BIN} setup.ini io_deisa.yml

simu_pid=$!
wait $simu_pid


