mkdir build bin
cd build
FC=gfortran cmake ..
make
cd ..
cp data/*.bsp bin/
<<<<<<< HEAD
cp data/naif0008.tls bin/
cp data/de414.bsp bin/
cd bin/ && ./ast_optim
=======
cp data/*.tls bin/
cd bin/ && clear && ./ast_optim
>>>>>>> 5867e6ab47fc84e2af5dcdc7f8ecf4773a9a7c27

