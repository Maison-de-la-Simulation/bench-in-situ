NDASKWORKERS=2


SCHEFILE="scheduler.json"

run_dask () {
  dask scheduler --dashboard-address localhost:8787 --scheduler-file=$SCHEFILE 1>scheduler.o 2>scheduler.e &
  while ! [ -f $SCHEFILE ]; do
      sleep 3
  done
  sync
  dask worker --dashboard-address localhost:8789 --nworkers $NDASKWORKERS --nthreads 1 --local-directory /tmp --scheduler-file=$SCHEFILE 1>worker.o 2>worker.e &
  sleep 3
}

# since we source this when calling scripts/nix/launch_*.sh, the relative path is relative to where 
# we call the script from, i.e the root directory.
SETUPINI="./simulation/setup.ini"
get_value() {
    local key=$1
    local file=$SETUPINI
    grep -E "^$key=" "$file" | cut -d'=' -f2 | tr -d ' '
}

# Read values from the ini file
mx=$(get_value "mx")
my=$(get_value "my")
mz=$(get_value "mz")

# remember to change simulation/setup.ini mx,my,mz to a combination so that mx * my * mz = NMPI
# Compute NMPI as the product of mx, my, and mz
NMPI=$((mx * my * mz))

echo "Parsed values from $SETUPINI: mx=$mx, my=$my, mz=$mz"
echo "Computed NMPI=$NMPI"
echo "NDASKWORKERS=$NDASKWORKERS"
echo "SCHEFILE=$SCHEFILE"

