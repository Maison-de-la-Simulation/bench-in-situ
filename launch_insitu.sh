SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
WORKING_DIR=${SCRIPT_DIR}/working_dir # TODO: Add datetime ?
SCHEFILE=${WORKING_DIR}/scheduler.json

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

dask-scheduler --scheduler-file=$SCHEFILE &
while ! [ -f $SCHEFILE ]; do
	sleep 3
done

dask-worker --local-directory /tmp --scheduler-file=$SCHEFILE &

python3 ../in-situ/fft.py

cd --
