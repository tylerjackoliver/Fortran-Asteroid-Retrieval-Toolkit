mkdir build bin
cd build
FC=gfortran cmake ..
make
cd ..
cd bin/ &&  ./pre_filter

