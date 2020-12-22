mkdir -p build bin
(
cd build
FC=mpif90 cmake ..
make
)
