#!/bin/bash


SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

source ${SCRIPT_DIR}/env.sh
source ${PDI_INSTALL_DIR}/share/pdi/env.sh

set -xeu

if [ -n $GUIX_ENVIRONMENT ]; then
  export LD_PRELOAD=$GUIX_ENVIRONMENT/lib/libz.so
fi

print_env

type pdirun
type python3

cd ${WORKING_DIR}

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

  source ${PYTHON_ENV}/bin/activate
  if [ -z "${VIRTUAL_ENV}" ]; then
    echo "Could not activate python environment !"
    exit 1
  fi

  run_dask
  #TODO: move client script to env/local folder
  python3 ${SCRIPT_DIR}/../../in-situ/fft.py $SCHEFILE &

elif [ "$1" = "new" ]; then
  echo "#####################"
  echo "# Running NEW Deisa #"
  echo "#####################"
  echo "Setting up Dask environement"

  PYTHON_ENV="${PYTHON_ENV}_new"
  source ${PYTHON_ENV}/bin/activate
  if [ -z "${VIRTUAL_ENV}" ]; then
    echo "Could not activate python environment !"
    exit 1
  fi

  run_dask
  #TODO: move client script to env/local folder
  python3 ${SCRIPT_DIR}/../../in-situ/bench_deisa.py ${DASK_NB_WORKERS} $SCHEFILE &
else
  echo "Unknown option. accepted values are: old, new"
  exit 1
fi


cp ${SCRIPT_DIR}/deisa/io.yml ${BUILD_DIR}
cp ${SCRIPT_DIR}/deisa/setup.ini ${BUILD_DIR}
cp ${SCHEFILE} ${BUILD_DIR}

cd ${BUILD_DIR}

echo "Running simulation"
#perf record mpirun -np 2 ./main setup.ini io_deisa.yml
mpirun -np ${MPI_NB_PROCS} ${SIMULATION_BIN} setup.ini io.yml

simu_pid=$!
wait $simu_pid

deactivate
pkill dask

