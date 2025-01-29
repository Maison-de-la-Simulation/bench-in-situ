#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env


run_python_import_tests () {
  echo "##############################"
  echo "# Testing Python environment #"
  echo "##############################"

  declare -a TO_IMPORT=("deisa" "h5py" "numpy" "dask" "distributed")
  
  for i in "${TO_IMPORT[@]}"
  do
    echo -n "import $i "
    python3 -c "import $i"
    retVal=$?
    if [ $retVal -ne 0 ]; then
      echo -e "${RED}Error${ENDCOLOR}"
      exit 1
    else
      echo -e "${GREEN}OK${ENDCOLOR}"
    fi
  done
}


cd ${WORKING_DIR}


echo "#######################"
echo "# Old Deisa"
echo "#######################"

# setup python environment
python3 -m venv ${PYTHON_ENV}

# activate python environment
source ${PYTHON_ENV}/bin/activate
pip install --upgrade pip

DEPS_FILE="deisa_deps_py${PYTHON_VERSION}.tar.gz"
if [ ! -f "$DEPS_FILE" ]; then
    echo "File not found: $DEPS_FILE in ${WORKING_DIR}"
    exit 1
fi

# untar and install dependencies
mkdir tmp_untar
tar xvf ${DEPS_FILE} -C tmp_untar
pip install tmp_untar/*.whl
rm -rf tmp_untar

# install Deisa
pip install --no-index --no-build-isolation --no-deps ../lib/deisa

run_python_import_tests

deactivate


echo "#######################"
echo "# New Deisa"
echo "#######################"
PYTHON_ENV="${PYTHON_ENV}_new"

# setup python environment
python3 -m venv ${PYTHON_ENV}

# activate python environment
source ${PYTHON_ENV}/bin/activate
pip install --upgrade pip

DEPS_FILE="new_deisa_deps_py${PYTHON_VERSION}.tar.gz"
if [ ! -f "$DEPS_FILE" ]; then
    echo "File not found: $DEPS_FILE in ${WORKING_DIR}"
    exit 1
fi

# untar and install dependencies
mkdir tmp_untar
tar xvf ${DEPS_FILE} -C tmp_untar
pip install tmp_untar/*.whl
rm -rf tmp_untar

# install Deisa
pip install --no-index --no-build-isolation --no-deps ../lib/new_deisa

run_python_import_tests

deactivate

cd --
