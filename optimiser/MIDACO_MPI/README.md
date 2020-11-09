# Parallel Mixed-Integer Distributed Ant Colony Optimiser


Uses [MIDACO](https://midaco-solver.com) to perform a global optimisation on the asteroid retrieval trajectories in invariant manifolds of the CR3BP.

The optimisation problem is defined by four continuous variables:

 * t<sub>0</sub>: Epoch of the transfer (Ephemeris seconds)
 * t<sub>t</sub>: Transfer time (s)
 * t<sub>end</sub>: Backwards integration time from the target plane (interpolated; dimensionless)
 * n<sub>mnfd</sub>: Discretisation along the target periodic orbit (interpolated; dimensionless)
 * K: Discretisation along the family of periodic orbits (interpolated; dimensionless)

This parallel version contains routines required optimise simultaneously on `m` cores.

This program requires support for the MPI parallel programming paradigm, from a minimum of version 3. The program has been tested using MPICH, OpenMPI and Intel MPI
on GNU and Intel compilers.

Any Fortran compilers must have support for the Fortran 2003 binding (in particular, `ISO_C_BINDING`.)

### Build Instructions

To build this application, alter the compiler used by changing the `FC=` environment variable to the command which invokes the correct MPI wrapper to your compiler of choice. For example, `FC=mpifort`, `FC=mpif90`, `FC=mpiifort`.

Once done, re-comment or de-comment the compiler options in `CMakeLists.txt` to match whether the GNU or Intel compiler is in use.

Running `build.sh` should automatically create the required folders, generate the makefiles and compile and link the program.

### Dependencies

* [jacobwilliams/bspline-fortran](https://github.com/jacobwilliams/bspline-fortran/) -- used for spline interpolation of the target conditions

* [jacobwilliams/Fortran-Astrodynamics-Toolkit](https://github.com/jacobwilliams/bspline-fortran/) -- used for general Lambert transfers 

* [NAIF Spice Toolkit](<https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html>) -- used for determining ephemeris data

* `CMAKE` & `make` -- used for building all of the source codes

* The codes written by `jacobwilliams` also depend on [FoBiS.py](https://github.com/szaghi/FoBiS), an automated build system designed specifically for Fortran, and written in Python.

Dependencies should be compiled for your machine, with the library files placed in a folder called `lib/` in each of the program directories (the Spice Toolkit library will need to be renamed from `spicelib.a` to `libspice.a` to be coherent with CMAKE conventions for library naming.) Any generated module files should be placed in a folder called `mod/` in each of the program directories.

### Usage instructions

Create a folder called data/ in the root directory of the optimiser. In it, the following files should be placed:

  * A valid planetary ephemeris file.
  * A valid leapseconds file.
  * Ephemeris files for any of the desired candidates; the get_ephemeris script found in this repository may be used to download the ephemeris files.
  * Data files for the manifold conditions; csv files with size n x 7, containing the backwards integration time in the synodic frame in column one, and the state vector in the remaining columns, where the units are in the synodic frame of the CR3BP.
  * Original orbit dataset: csv file with size n x 6, containing the state vector of initial conditions, where the units are in the synodic frame of the CR3BP.
  
The file `problem_parameters.f90` can then be adjusted to your particular problem:

```Fortran
module problem_parameters
    
     

    implicit none

    ! ////////////////////////////////////////////////////////////
    !
    ! Contains variables that control the program runtime
    !
    ! Inputs
    ! ~~~~~~
    ! targ_can: SPK ID of the target candidate
    ! datafile: Location of target dataset file
    ! dataset: Allocatable double-precision array that contains the target section
    ! is_loaded: Whether the dataset has been loaded into memory
    ! minimum_transfer_time: Minimum duration of Lambert arc (optimiser constraint)
    ! maximum_transfer_time: Maximum duration of Lambert arc (optimiser constraint)
    !
    ! ///////////////////////////////////////////////////////////

    ! //////////////////////////////
    ! 
    ! PROBLEM PARAMETERS - CHANGE ME
    !
    character(7)                    :: targ_can                                                             ! Needs to change to fit the target *currently* being studied
    character(*),  parameter        :: datafile = &
        '../path/to/states/at/the/pi/8/section/'

    character(*),  parameter        :: original_orbit_data = &
        '../path/to/original/orbit/conditions/'

    character(*),  parameter        :: target_file = '../data/spk_ids_to_optimise'

    character(*),  parameter        :: ephemeris_prefix=&
        '../path/to/ephemeris/files/'

    double precision, parameter     :: minimum_transfer_time = 0.d0* 86400                                  ! Seconds
    double precision, parameter     :: maximum_transfer_time = 1500.d0 * 86400                              ! Seconds
    integer,       allocatable      :: targ_can_array(:)                                                    ! Placeholder for now
    integer,       parameter        :: t_end_disc = 1000                                                    ! Number of backwards-time integrations
    integer,       parameter        :: n_mnfd_disc = 360                                                    ! Number of in-plane discretisations
    integer,       parameter        :: max_time = 10 * 60 ! 5 hours
    ! //////////////////////////////

    real*8, POINTER                 :: dataset(:,:,:) => NULL()                                             ! t_end x n_mnfd x J x 6
    real*8, POINTER                 :: perturbed_conds_dataset(:,:,:) => NULL()                             ! n_mnfd x J x 6

    double precision                :: transfer_time
    double precision                :: transfer_time_array(1500)
    double precision                :: pareto_front_minimum

end module problem_parameters
```

### Tested compilers and platforms

MacOS, gfortran-8, gcc-8, ifort, intel-mpi, OpenMPI, MPICH, Python 3.6.

Ubuntu, gfortran-7, gcc-7, ifort, intel-mpi, OpenMPI, MPICH, Python 3.4.
Red Hat, gfortran-7, gcc-7, ifort, intel-mpi, OpenMPI, MPICH, Python 3.4.
