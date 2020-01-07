mkdir build bin
cd build
FC=gfortran cmake ..
make
cd ..
cp data/*.bsp bin/
cp data/naif0008.tls bin/
cp data/de414.bsp bin/
cd bin/ && clear && echo "3550232" && ./ast_optim

