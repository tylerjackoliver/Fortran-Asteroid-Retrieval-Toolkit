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
