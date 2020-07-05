mkdir build bin
cd build
FC=mpif90 cmake ..
make
cd ../bin
mpirun -np 4 ./ast_optim