#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env


if [ "$1" = "old" ]; then
  echo "#####################"
  echo "# Running OLD Deisa #"
  echo "#####################"

  source ${WORKING_DIR}/${PYTHON_ENV}/bin/activate
  source ${PDI_INSTALL_DIR}/share/pdi/env.sh

  #export PYTHONPATH=${WORKING_DIR}/${PYTHON_ENV}/lib/python3.11/site-packages:${PYTHONPATH}
  #export LD_LIBRARY_PATH=${PDI_INSTALL_DIR}/lib:${WORKING_DIR}/${PYTHON_ENV}/lib:${LD_LIBRARY_PATH}
  #export LD_LIBRARY_PATH=${WORKING_DIR}/${PYTHON_ENV}/lib:${LD_LIBRARY_PATH}
  #export LD_PRELOAD=${PDI_INSTALL_DIR}/lib64/libyaml.so # note: only for mac
elif [ "$1" = "new" ]; then
  echo "#####################"
  echo "# Running NEW Deisa #"
  echo "#####################"
  
  source ${WORKING_DIR}/${PYTHON_ENV}_new/bin/activate
  source ${PDI_INSTALL_DIR}/share/pdi/env.sh

  #export PYTHONPATH=${WORKING_DIR}/${PYTHON_ENV}/lib/python3.11/site-packages:${PYTHONPATH}
  #export LD_LIBRARY_PATH=${PDI_INSTALL_DIR}/lib:${WORKING_DIR}/${PYTHON_ENV}/lib:${LD_LIBRARY_PATH}
  #export LD_LIBRARY_PATH=${WORKING_DIR}/${PYTHON_ENV}/lib:${LD_LIBRARY_PATH}
  #export LD_PRELOAD=${PDI_INSTALL_DIR}/lib64/libyaml.so # note: only for mac
else
  echo "Unknown option. accepted values are: old, new"
  exit 1
fi





echo "PYTHONPATH=${PYTHONPATH}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
echo "LD_PRELOAD=${LD_PRELOAD}"

cp ../simulation/io_chkpt.yml ../simulation/io_deisa.yml ${BUILD_DIR}
cp ../simulation/setup.ini ${BUILD_DIR}
cp ${SCHEFILE} ${BUILD_DIR}

cd ${BUILD_DIR}

#perf record mpirun -np 2 ./main setup.ini io_deisa.yml
mpirun -np 4 ./main setup.ini io_deisa.yml

#pdirun mpirun -np 1 ${BUILD_DIR}/main ../setup.ini ../io_deisa.yml
cd --
