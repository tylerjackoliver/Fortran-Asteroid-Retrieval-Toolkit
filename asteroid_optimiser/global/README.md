# Global optimiser

Contains routines to perform the global optimiser on the Lambert transfers.

## Program set-up

Create a folder called data/ in the root directory of the optimiser. In it, the following files should be placed:

  * de414.bsp
  * naif0008.tls
  * Ephemeris files for any of the desired candidates.
  * Data files for the manifold conditions; csv files with size n x 6, where the units are km, km/s in the inertial frame.
  
A proper interface to this function is yet to be developed. As such, the filename for the datafile must be manually edited in ```main/src.f90```, as well as all references to the target candidate. A COMMON block to handle this data will be implemented in the near future.

## Building

These subroutines use CMake to build the necessary softwares.

This program relies on the proper libraries for the SPICE toolkit in the proper compiler and platform version being present in lib/ (and renamed to be libspice.a, since the default name of spicelib.a is incompatible with default CMAKE library names.)

The [Fortran Astrodynamics Toolkit](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit) is required for the Lambert solver. FAT itself requires [FoBiS.py](https://github.com/szaghi/FoBiS), a Fortran-specific proprietary automatic build system; for all its flaws, it's the only way to get the FAT working. The module files for FAT should be placed in mod/.

To build, simply use:

```
sh build.sh
```

at the terminal. Similar scripts exist to clean the directory, and rebuild (clean + build).

## Tested compilers and platforms

MacOS, gfortran-8, gcc-8, Python 3.6.
Ubuntu, gfortran-7, gcc-7, Python 3.4.
