#!/bin/bash

source env.sh
print_env

cd ${WORKING_DIR}

dask scheduler --scheduler-file=$SCHEFILE &
sleep 3
sync

dask worker --local-directory /tmp --scheduler-file=$SCHEFILE &

python3 ../in-situ/fft.py $SCHEFILE

cd --
