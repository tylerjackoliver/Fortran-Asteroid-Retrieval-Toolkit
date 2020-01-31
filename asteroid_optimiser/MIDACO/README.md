# Mixed-Integer Distributed Ant Colony Optimiser

Uses [MIDACO](https://midaco-solver.com) to perform a global optimisation on the asteroid retrieval trajectories in invariant manifolds of the CR3BP.

The optimisation problem is defined by three continuous variables:

 * t<sub>0</sub>: Epoch of the transfer (Ephemeris seconds)
 * t<sub>t</sub>: Transfer time (s)
 * t<sub>end</sub>: Backwards integration time from the target plane (interpolated; dimensionless)
 * n<sub>mnfd</sub>: Discretisation along the target periodic orbit (interpolated; dimensionless)

## Program set-up

Create a folder called data/ in the root directory of the optimiser. In it, the following files should be placed:

  * de414.bsp
  * naif0008.tls
  * Ephemeris files for any of the desired candidates.
  * Data files for the manifold conditions; csv files with size n x 7, containing the backwards integration time in the synodic frame in column one, and the state vector in the remaining columns, where the units are km, km/s in the inertial frame.
  * Original orbit dataset: csv file with size n x 7, containing the backwards integration time in the synodic frame in column one, and the state vector in the remaining columsns, where the units are km, km/s in the inertial frame.
  
The file `problem_parameters.f90` can then be adjusted to your particular problem:

```Fortran

    ! ///////////////////////////////////////////////////////////
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
    character(*),  parameter :: targ_can = '3390109'
    character(*),  parameter :: data_file = '../data/2020-01-12_L2PlanarBackCondsGlobal.csv'
    character(*),  parameter :: original_orbit_data = '../data/original_orbits.csv'
    real(kind=dp), parameter :: minimum_transfer_time = 1.d0* 86400                          ! Seconds
    real(kind=dp), parameter :: maximum_transfer_time = 1600.d0 * 86400                      ! Seconds
    !
    ! //////////////////////////////

```

If further tweaking is required (i.e. the input data-file is discretized differently), then `variable_initialisation.f90` also needs to be edited, particularly the discretisation measures:

```Fortran

    integer, parameter      :: t_end_disc = 100                                             ! Parameterisation in backwards time
    integer, parameter      :: n_mnfd_disc = 360                                            ! Parameterisation around the orbit

```

## Building

To build, simply use:

```
sh build.sh
```

at the terminal. Similar scripts exist to clean the directory, and rebuild (clean + build).

## Tested compilers and platforms

MacOS, gfortran-8, gcc-8, Python 3.6.

Ubuntu, gfortran-7, gcc-7, Python 3.4.
