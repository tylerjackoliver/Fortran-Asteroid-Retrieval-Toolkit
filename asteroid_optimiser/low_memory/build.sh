mkdir build bin
cd build
FC=mpif90 cmake ..
make
cd ..
cp data/*.bsp bin/
cp data/*.tls bin/
cd bin/ && mpirun -np 4 ./ast_optim

