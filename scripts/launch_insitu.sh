#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env


run_dask () {
  dask scheduler --scheduler-file=$SCHEFILE &
  sleep 3
  sync
  dask worker --local-directory /tmp --scheduler-file=$SCHEFILE &
}



cd ${WORKING_DIR}

if [ "$1" = "old" ]; then
  echo "#####################"
  echo "# Running OLD Deisa #"
  echo "#####################"

  source ${PYTHON_ENV}/bin/activate
  run_dask
  python3 -O ../in-situ/fft.py $SCHEFILE
elif [ "$1" = "new" ]; then
  echo "#####################"
  echo "# Running NEW Deisa #"
  echo "#####################"
  
  source ${PYTHON_ENV}_new/bin/activate
  run_dask
  python3 -O ../in-situ/bench_deisa.py $SCHEFILE
else
  echo "Unknown option. accepted values are: old, new"
  exit 1
fi


cd --
