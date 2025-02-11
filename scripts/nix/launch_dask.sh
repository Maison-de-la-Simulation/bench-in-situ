#!/bin/bash

set -xeu

OPWD=$PWD

# set NDASKWORKERS, NMPI, SCHEFILE, and run_dask
source $PWD/scripts/nix/env.sh

EXPERIMENT_DIR=$PWD/experiment

if [ -d "$EXPERIMENT_DIR" ]; then
    echo "Running new experiment... Deleting old ones."
    rm -rf $EXPERIMENT_DIR
    mkdir -p $EXPERIMENT_DIR
else
    echo "Creating new experimemnt directory..."
    mkdir -p $EXPERIMENT_DIR
fi

cd $EXPERIMENT_DIR

run_dask
cd $OPWD
set +xeu
