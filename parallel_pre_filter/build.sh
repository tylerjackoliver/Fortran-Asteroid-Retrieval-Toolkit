mkdir build bin
cd build
FC=mpif90 cmake ..
make
cd ..
cd bin/ && mpirun -np 4 ./pre_filter
