mkdir -p build bin
(
cd build
FC=mpif90 cmake ..
make
)
cp data/*.bsp bin/
cp data/*.tls bin/
(
cd bin/
mpirun -np 2 ./ast_optim
)