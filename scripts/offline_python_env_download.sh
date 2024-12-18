#!/bin/bash

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd) && cd ${SCRIPT_DIR}
source env.sh
print_env

py_versions=( "3.9" )
for py_version in "${py_versions[@]}"; do
    echo "Working on python version: ${py_version}"
    docker run -it --rm -w /root/workingdir -v ${SCRIPT_DIR}/..:/root/workingdir python:${py_version} bash -c "pip install --upgrade pip && \
    pip download lib/deisa -d ${PYTHON_WORKING_DIR} && \
    pip download wheel versioneer zarr cytools h5py -d ${PYTHON_WORKING_DIR} && \
    cd ${PYTHON_WORKING_DIR} && pip wheel ${PYTHON_WORKING_DIR}/distributed-2021.11.2+1398.g9b7ce185.zip && cd - && \
    rm ${PYTHON_WORKING_DIR}/distributed-2011.11.2+1398.g9b7ce185.zip && \
    tar -cvzf output/deisa_deps_py${py_version//.}.tar.gz -C ${PYTHON_WORKING_DIR} . && \
    rm -rf ${PYTHON_WORKING_DIR}"
done

