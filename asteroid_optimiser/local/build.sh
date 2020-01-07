mkdir build bin
cd build
FC=gfortran cmake ..
make
cd ..
cp data/*.bsp bin/
cp data/naif0008.tls bin/
cp data/de414.bsp bin/
cp data/2019-11-20_L2PlanarBackCondsGlobal.csv bin/
cd bin/ && echo "3550232 local 332" && ./ast_optim

