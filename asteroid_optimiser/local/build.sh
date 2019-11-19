mkdir build bin
cd build
export CC=/usr/local/bin/gcc-8
export CXX=/usr/local/bin/g++-8
FC=gfortran cmake ..
make
cd ..
cp data/3390109.bsp bin/
cp data/naif0008.tls bin/
cp data/de414.bsp bin/
cp data/test_conditions_500000.dat bin/
cd bin/ && ./ast_optim

