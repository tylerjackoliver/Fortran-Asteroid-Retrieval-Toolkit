mkdir build bin
cd build
FC=gfortran cmake ..
make
cd ..
cp data/*.bsp bin/
cp data/*.tls bin/
cd bin/ && ./ast_optim
