mkdir build bin
cd build
export CC=/usr/local/bin/gcc-8
export CXX=/usr/local/bin/g++-8
FC=gfortran cmake ..
make
cd ..
cp data/*.bsp bin/
cp data/*.tls bin/
cd bin/ && ./ast_optim
