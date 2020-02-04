module problem_parameters
    
    use precision_kinds

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
    character(*),  parameter :: targ_can = '3390109'
    character(*),  parameter :: datafile = '../data/2020-01-28_L2PlanarBackCondsGlobal.csv'
    character(*),  parameter :: original_orbit_data = '../data/2020-02-01_L2PlanarPerturbedConds.csv'
    real(kind=dp), parameter :: minimum_transfer_time = 1.d0* 86400                        ! Seconds
    real(kind=dp), parameter :: maximum_transfer_time = 1600.d0 * 86400                      ! Seconds
    !
    ! //////////////////////////////

    real*8,     allocatable :: dataset(:,:,:,:)                                             ! t_end x n_mnfd x J x 6
    real*8,     allocatable :: perturbed_conds_dataset(:,:,:)                               ! n_mnfd x J x 6

end module problem_parameters
