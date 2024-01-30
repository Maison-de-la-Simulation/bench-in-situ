#git clone --recurse-submodules https://github.com/Maison-de-la-Simulation/bench-in-situ.git
#cd bench-in-situ

rm -r bench-in-situ
git clone https://github.com/Maison-de-la-Simulation/bench-in-situ.git
cd bench-in-situ

mkdir -p vendors && cd vendors

cd ..
cd lib
wget https://github.com/kokkos/kokkos/archive/refs/tags/4.2.00.tar.gz
gunzip 4.2.00.tar.gz
tar -xf 4.2.00.tar

wget https://github.com/benhoyt/inih/archive/refs/tags/r58.tar.gz
gunzip r58.tar.gz
tar -xf r58.tar
