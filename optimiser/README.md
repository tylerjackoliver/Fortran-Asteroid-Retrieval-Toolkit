# Asteroid optimisers

## Current status

The main routine is now MIDACO/.

## Overview

* MIDACO/ uses the MIDACO-solver with interpolation to attempt to locate the minima of the objective function; the current main optimiser, it has been tested and is designed to automatically return the Pareto front and interface with MATLAB validation routines.
* global/ implements a Derivative-Free Boender-Timmer-Rinnoy Kan algorithm to attempt to find the global minimum. 
* local/ implements the BOBYQA derivative-free local optimiser. Not often used; rarely makes improvements on the global solver.
* ephemeris/ implements the global solver but by using ephemerides at all points of the optimisation. Not refactored.
* global_test/ is a `playground' for testing various things in the optimisation.

## Dependencies

Specific dependencies are listed in the README.md of every subfolder; however, the majority of the programs rely on the modules below:

* [jacobwilliams/bspline-fortran](https://github.com/jacobwilliams/bspline-fortran/) -- used for spline interpolation of the target conditions

* [jacobwilliams/Fortran-Astrodynamics-Toolkit](https://github.com/jacobwilliams/bspline-fortran/) -- used for general Lambert transfers 

* [NAIF Spice Toolkit](<https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html>) -- used for determining ephemeris data

* `CMAKE` & `make` -- used for building all of the source codes

* The codes written by `jacobwilliams` also depend on [FoBiS.py](https://github.com/szaghi/FoBiS), an automated build system designed specifically for Fortran, and written in Python.

* A good Fortran compiler, and an MPI wrapper to it; these codes are written and tested using `gfortran-8` and `ifort`.

Dependencies should be compiled for your machine, with the library files placed in a folder called `lib/` in each of the program directories (the Spice Toolkit library will need to be renamed from `spicelib.a` to `libspice.a` to be coherent with CMAKE conventions for library naming.) Any generated module files should be placed in a folder called `mod/` in each of the program directories.
