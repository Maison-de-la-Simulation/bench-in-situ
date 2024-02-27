
1) Copy the whole repository to your home directory on Ruche (i.e.: /gpfs/users/<USERNAME>).
2) Copy the `deisa_deps_pyXX.tar.gz` archive into the `working_dir` folder (as defined in `env.sh`)
3) From the copied folder, run the `build_all.sh` script: `./scipts/build_all.sh`. This will setup the python environment for Deisa, build PDI and the simulation.
4) Submit Slurm tasks using the slurm script: `sbatch envs/ruche/launcher.sh`



