mkdir build bin
cd build
FC=gfortran cmake ..
make
cd ..
cp data/3435539.bsp bin/
cp data/3390109.bsp bin/
cp data/naif0008.tls bin/
cp data/de414.bsp bin/
cd bin/ && ./ast_optim

