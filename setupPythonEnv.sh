#python -m venv ./PyEnv
source PyEnv/bin/activate
cd pyLibs

tar -xzvf PyYAML-6.0.1.tar.gz
tar -xzvf numpy-1.26.3.tar.gz
tar -xzvf dask-2023.11.0.tar.gz
tar -xzvf Cython-3.0.8.tar.gz
tar -xzvf setuptools-69.0.3.tar.gz
tar -xzvf wheel-0.42.0.tar.gz
tar -xzvf flit_core-3.9.0.tar.gz

rm -r *.tar.gz

cd flit_core-3.9.0
pip install .

cd ..
cd wheel-0.42.0
pip install .

cd ..
cd setuptools-69.0.3
pip install .

cd ..
cd Cython-3.0.8
pip install .

cd distributed
pip install .

cd ..
cd numpy-1.26.3
pip install .

cd ..
cd PyYAML-6.0.1
pip install .

cd ..
cd dask-2023.11.0
pip install .


