mkdir build bin
cd build
FC=gfortran cmake ..
make
cd ..
cp data/*.bsp bin/
cp data/naif0008.tls bin/
cp data/de414.bsp bin/
cd bin/ && ./ast_optim

