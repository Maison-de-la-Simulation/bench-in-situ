# WARNING

At the moment the deisa is being tested and fixed. After running 
`git submodule init && git submodule update` navigate to `./lib/new_deisa` and 
run `git switch fix_graceful_shutdown`. This will put you in the correct branch with 
the fix. When the fix is merged in the main branch, I will remove this. For now I leave it at 
the top so its the first thing I see when reading this (to remind me to deleted it).


# Launch Scripts When Using Nix

This directory contains the scripts needed to run the benchmarks when using Nix to set up the 
environment (only for the new version of deisa). The logic is simpler than the other ones and 
relies on Nix completely setting up the environment. **All scripts should be launched from the 
root path of the project**.

## Setup Nix Environment

Before running the scripts, run `git submodule init` and `git submodule update` to clone the 
repositories of the projects. 

Be sure to point the environment variable `PYTHONPATH` to the root of the Deisa project.  
Also, be sure to activate the nix environment to get all other dependencies.
```bash
# If using git submodule it should be ROOT_OF_PROJECT/libs/new_deisa
$ export PYTHONPATH=ROOT_OF_DEISA 
$ nix-shell ROOT_OF_PROJECT/envs/nix/shell.nix
```
(Optional) Test that python finds deisa by running `import deisa` from within the interpreter.
(Optional) A convenient `.envrc` is located at the root of the project for people using `direnv`.

## Building The Simulation

To build the simulation run:
```bash
$ bash ROOT_OF_PROJECT/scripts/nix/build_simulation.sh
```
## Setting Up Parameters

The main entry point for all parameters is `scripts/envs/nix/env.sh` and `./simulation/setup.ini`:
- `scripts/envs/nix/env.sh`: reads `setup.ini` (look at the next bullet point) to set `NMPI`, sets 
the number of dask workers `NDASKWORKERS`, and the name of the scheduler file `SCHEFILE`.
- `./simulation/setup.ini`: controls various simulation-related parameters. The number of MPI 
processes used (`NMPI`) to run the simulation is be the product of `mx`, `my` and `mz` which are 
set in `setup.ini` under `[mesh]`.

## Launching Analytics And Simulation

To launch the benchmark, you can follow two paths: 
```bash
# launch everything from a single script - great if you want to run 
# everything once.
$ bash ROOT_OF_PROJECT/scripts/nix/launch_all.sh

# or 

# setup the cluster first and then launch the analytics + simulation 
# this is great if you want to run the benchmark many times as you modify it since 
# it reuses the same dask cluster
$ bash ROOT_OF_PROJECT/scripts/nix/launch_dask.sh
$ bash ROOT_OF_PROJECT/scripts/nix/launch_analytics_sim.sh
```
The results of the experiment run will be placed in a newly created `experiment` directory at the 
root of the project. 

**This directory is eliminated each time you run `launch_dask.sh` or `launch_all.sh`**. Please modify 
the scripts accordingly if you want this behavior to change to preserve all runs.

# Prometheus Monitoring (Optional)

The nix shell also setups the necessary environment for prometheus monitoring. This is useful 
since it provides some additional diagnostics. 

To run a prometheus instance, run `launch_dask.sh` and then run 
`prometheus --config.file=prometheus.yaml` where `prometheus.yaml` is copied from the `envs/nix` 
directory (or simply point to it). At the moment, it only tracks the metrics exposed by the dask 
scheduler dashboard (`localhost:8787`) and the dask worker dashboard (`localhost:8789`) - both 
expose prometheus analytics at `/metrics` endpoint. To view the prometheus metrics, simply go 
to `localhost:9090`.

