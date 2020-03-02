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
    character(7)                :: targ_can                                                             ! Needs to change to fit the target *currently* being studied
    character(*),  parameter    :: datafile = '../data/2020-02-14_L1PlanarBackCondsSynodicDense.csv'
    character(*),  parameter    :: original_orbit_data = '../data/2020-02-01_L2PlanarPerturbedConds.csv'
    character(*),  parameter    :: target_file = '../data/targets_to_consider'
    real(kind=dp), parameter    :: minimum_transfer_time = 1.d0* 86400                                  ! Seconds
    real(kind=dp), parameter    :: maximum_transfer_time = 1600.d0 * 86400                              ! Seconds
    integer,       allocatable  :: targ_can_array(:)                                                    ! Placeholder for now
    integer,       parameter    :: t_end_disc = 1000                                                    ! Number of backwards-time integrations
    integer,       parameter    :: n_mnfd_disc = 360                                                    ! Number of in-plane discretisations
    !
    ! //////////////////////////////

    real*8, POINTER             :: dataset(:,:,:,:) => NULL()                                           ! t_end x n_mnfd x J x 6
    real*8, POINTER             :: perturbed_conds_dataset(:,:,:) => NULL()                             ! n_mnfd x J x 6

end module problem_parameters
