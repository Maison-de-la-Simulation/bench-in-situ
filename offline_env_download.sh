#!/bin/bash

source env.sh
print_env

mkdir -p ${WORKING_DIR}
cd ${WORKING_DIR}

# Create conda python 3.9
# 1) create conda env with python 3.9: conda create -y --name py39 python=3.9.10
conda create -y --name py39 python=3.9
# 2) activate environment
conda init          #TODO: can this be removed ?
conda activate py39 #TODO: use `conda run` instead of `conda activate`
# 3) pip download and wheel
pip install ../lib/deisa            #TODO: is this needed ?
pip download ../lib/deisa/ -d deps/ #TODO: deps should be a a variable in env.sh
pip download wheel versioneer -d deps/

cd deps
pip wheel *.zip
rm *.zip
cd ..

# 4) tar gz -> deisa_deps_py39.tar.gz
tar -czvf deisa_deps_py39.tar.gz deps/*
# 5) cleanup: deactivate conda env, rm conda env, downloaded files
conda deactivate
conda remove -y --name py39 --all
rm -r deps/*

# Create conda python 3.10

# 1) create conda env with python 3.10: conda create -y --name py39 python=3.9.10
conda create -y --name py10 python=3.10
# 2) activate environment
conda init
conda activate py310
# 3) pip download and wheel
pip install ../lib/deisa
pip download ../lib/deisa/ -d deps/
pip download wheel versioneer -d deps/

cd deps
pip wheel *.zip
rm *.zip
cd ..

# 4) tar gz -> deisa_deps_py39.tar.gz
tar -czvf deisa_deps_py310.tar.gz deps/*
# 5) cleanup: deactivate conda env, rm conda env, downloaded files
conda deactivate
conda remove -y --name py310 --all
rm -r deps/*

# Create conda python 3.11
# 1) create conda env with python 3.11: conda create -y --name py39 python=3.9.10
conda create -y --name py311 python=3.11
# 2) activate environment
conda init
conda activate py311
# 3) pip download and wheel
pip install ../lib/deisa
pip download ../lib/deisa/ -d deps/
pip download wheel versioneer -d deps/

cd deps
pip wheel *.zip
rm *.zip
cd ..

# 4) tar gz -> deisa_deps_py39.tar.gz
tar -czvf deisa_deps_py311.tar.gz deps/*
# 5) cleanup: deactivate conda env, rm conda env, downloaded files
conda deactivate
conda remove -y --name py311 --all
rm -r deps/*

# setup python environment
#python -m venv ${PYTHON_ENV}

# activate python environment
#source ${PYTHON_ENV}/bin/activate

# install Deisa
#pip install ../lib/deisa

#pip download ../lib/deisa/ -d deps/
#pip download wheel versioneer -d deps/

#deactivate
rm -r ${PYTHON_ENV}
cd ..
