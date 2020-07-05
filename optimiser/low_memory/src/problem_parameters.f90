module problem_parameters
    
    implicit none

    ! ////////////////////////////////////////////////////////////
    !
    ! Contains variables that control the program runtime
    !
    ! ///////////////////////////////////////////////////////////

    character(7)                    :: targ_can                                                             ! Candidate we're studying
    
    character(*),  parameter        :: datafile = &
   '/media/jack/Files/asteroid_retrieval_datasets/Plane/planeConds.csv'                                     ! Target manifold dataset

    character(*),  parameter        :: original_orbit_data = &
    '/media/jack/Files/asteroid_retrieval_datasets/Perturbed/2020-03-03_L2VertLyapPerturbedConds.csv'       ! Initial perturbed orbit conditions

    character(*),  parameter        :: target_file = '../data/targets_to_consider'                          ! List of SPK IDs to optimise

    character(*),  parameter        :: ephemeris_prefix=&  
    '/media/jack/Files/ephemeris_files/ephemeris_files/'                                                    ! Location of relevant epehemeris files

    ! Optimisation settings

    double precision, parameter     :: minimum_transfer_time = 1.d0* 86400                                  ! Seconds
    double precision, parameter     :: maximum_transfer_time = 1400.d0 * 86400                              ! Seconds

    integer,          parameter     :: max_time = 1 * 60                                                    ! Amount of time to perform the global optimisation (s) (MIDACO)

    integer,          parameter     :: n_mnfd_disc = 360                                                    ! Number of in-plane discretisations

    real*8,           POINTER       :: dataset(:,:,:) => NULL()                                             ! n_mnfd x J x 7
    real*8,           POINTER       :: perturbed_conds_dataset(:,:,:) => NULL()                             ! n_mnfd x J x 6
    
    integer,       allocatable      :: targ_can_array(:)                                                    ! Placeholder for now
    
end module problem_parameters
