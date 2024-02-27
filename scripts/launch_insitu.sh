#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env

cd ${WORKING_DIR}

dask scheduler --scheduler-file=$SCHEFILE &
sleep 3
sync

dask worker --local-directory /tmp --scheduler-file=$SCHEFILE &

python3 -O ../in-situ/fft.py $SCHEFILE

cd --
