# Parallel Asteroid Pre-Filtering Tool

This program implements a shared-memory-parallelized version of the Fortran Asteroid Prefilter (FAP).

Uses [MIDACO](https://midaco-solver.com) to perform a pre-filtering scan on candidate asteroids from the Minor Bodies Database.

### Build Instructions

To build this application, alter the compiler used by changing the `FC=` environment variable to the command which invokes the correct MPI wrapper to your compiler of choice. For example, `FC=mpifort`, `FC=mpif90`, `FC=mpiifort`.

Once done, re-comment or de-comment the compiler options in `CMakeLists.txt` to match whether the GNU or Intel compiler is in use.

Running `build.sh` should automatically create the required folders, generate the makefiles and compile and link the program.

### Dependencies

* [NAIF Spice Toolkit](<https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html>) -- used for determining ephemeris data

* `CMAKE` & `make` -- used for building all of the source codes

* The codes written by `jacobwilliams` also depend on [FoBiS.py](https://github.com/szaghi/FoBiS), an automated build system designed specifically for Fortran, and written in Python.

Dependencies should be compiled for your machine, with the library files placed in a folder called `lib/` in each of the program directories (the Spice Toolkit library will need to be renamed from `spicelib.a` to `libspice.a` to be coherent with CMAKE conventions for library naming.) Any generated module files should be placed in a folder called `mod/` in each of the program directories.

### Usage instructions

Create a folder called data/ in the root directory of the optimiser. In it, the following files should be placed:

  * de414.bsp
  * naif0008.tls
  * Ephemeris files for any of the desired candidates.
  * Data files for the manifold conditions; csv files with size n x 7, containing the backwards integration time in the synodic frame in column one, and the state vector in the remaining columns, where the units are in the synodic frame of the CR3BP

The file `problem_parameters.f90` can then be adjusted to your particular problem.

### Tested compilers and platforms

MacOS, gfortran-8, gcc-8, ifort, intel-mpi, OpenMPI, MPICH, Python 3.6.

Ubuntu, gfortran-7, gcc-7, ifort, intel-mpi, OpenMPI, MPICH, Python 3.4.
Red Hat, gfortran-7, gcc-7, ifort, intel-mpi, OpenMPI, MPICH, Python 3.4.
