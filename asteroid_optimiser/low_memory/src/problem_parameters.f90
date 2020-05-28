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
    ! minimum_transfer_time: Minimum duration of Lambert arc (optimiser constraint)
    ! maximum_transfer_time: Maximum duration of Lambert arc (optimiser constraint)
    !
    ! ///////////////////////////////////////////////////////////

    character(7)                    :: targ_can 
    
    character(*),  parameter        :: datafile = &
   '/media/jack/Files/asteroid_retrieval_datasets/Plane/planeConds.csv'

    character(*),  parameter        :: original_orbit_data = &
    '/media/jack/Files/asteroid_retrieval_datasets/Perturbed/2020-03-03_L2VertLyapPerturbedConds.csv'

    character(*),  parameter        :: target_file = '../data/targets_to_consider'

    character(*),  parameter        :: ephemeris_prefix=&
        '/media/jack/Files/ephemeris_files/ephemeris_files/'

    ! Optimisation settings

    double precision, parameter     :: minimum_transfer_time = 1.d0* 86400                                  ! Seconds
    double precision, parameter     :: maximum_transfer_time = 1400.d0 * 86400                              ! Seconds

    integer,          parameter     :: max_time = 1 * 60                                                    ! Amount of time to perform the global optimisation (s)
    integer,          parameter     :: max_time_local = 1 * 60                                              ! Amount of time to perform the local optimisation * per solution * 

    integer,          parameter     :: t_end_disc = 1000                                                    ! Number of backwards-time integrations
    integer,          parameter     :: n_mnfd_disc = 360                                                    ! Number of in-plane discretisations

    real*8,           POINTER       :: dataset(:,:,:) => NULL()                                             ! n_mnfd x J x 6
    real*8,           POINTER       :: perturbed_conds_dataset(:,:,:) => NULL()                             ! n_mnfd x J x 6
    
    integer,       allocatable      :: targ_can_array(:)                                                    ! Placeholder for now
    
end module problem_parameters
