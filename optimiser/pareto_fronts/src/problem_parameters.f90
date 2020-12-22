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
    character(7)                        :: targ_can                                                             ! SPK ID of the target
    character(*),       parameter       :: datafile = '../path/to/datafile'                                     ! File with manifold reference states on +- pi/8 plane
    character(*),       parameter       :: target_file = '../data/targets_to_consider'                          ! Text file containing line-delimited SPK IDs of asteroids to study
    character(*),       parameter       :: ephemeris_prefix='../path/to/ephemeris/'                             ! Parent folder containing ephemeris for asteroids in target_file

    double precision,   parameter       :: minimum_transfer_time = 0.d0* 86400                                  ! Lower bound for transfer time, seconds
    double precision,   parameter       :: maximum_transfer_time = 1500.d0 * 86400                              ! Upper bound for transfer time, seconds
    integer,            parameter       :: n_mnfd_disc = 360                                                    ! Number of in-plane discretisations
    integer,            parameter       :: max_time = 10 * 60                                                   ! Maximum optimisation time, seconds
    ! //////////////////////////////

    real*8, POINTER                     :: dataset(:,:,:) => NULL()                                             ! Main reference state dataset

    integer,            allocatable     :: targ_can_array(:)                                                    ! All the targets to study

    double precision                    :: transfer_time                                                        ! Current transfer time
    double precision                    :: transfer_time_array(1500)                                            ! All transfer times
    double precision                    :: pareto_front_minimum                                                 ! Minimum value in the PF

end module problem_parameters
