#!/bin/bash
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
cd ${SCRIPT_DIR}

BASE_DIR=${SCRIPT_DIR}

RESULT_DIR=$1
RESULT_FILE="bench_output_"$1".txt"

echo $1 > ${BASE_DIR}/${RESULT_DIR}/${RESULT_FILE}
cat ${SCRIPT_DIR}/../../../lib/pdi/pdi/VERSION >> ${BASE_DIR}/${RESULT_DIR}/${RESULT_FILE}

cd ${BASE_DIR}/${RESULT_DIR} 
#grep -rw . -e "RESULT" >> ${BASE_DIR}/${RESULT_DIR}/${RESULT_FILE}
## To avoid warning
find . -name res*.out | xargs grep -e RESULT >>  ${BASE_DIR}/${RESULT_DIR}/${RESULT_FILE}

### create my "result.txt" ici

### case weak
cd ${BASE_DIR}/${RESULT_DIR}
python3 ../plotter_weak_scaling.py --input ${RESULT_FILE}

### case strong
#cd ${BASE_DIR}/${RESULT_DIR}
#python3 ../plotter_cube_sizes.py --input ${RESULT_FILE}


