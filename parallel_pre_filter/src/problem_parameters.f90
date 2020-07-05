module problem_parameters

    use precision_kinds
    use constants

    implicit none

    character(*), parameter                 :: target_dataset = '../data/your_dataset_here.csv'                         ! Input manifold conditions. num_targets x 7 (time, state - synodic)
    character(*), parameter                 :: spk_database   = '../data/your_targets_here'                             ! List of SPK designations to consider in the pre-filter
    character(*), parameter                 :: week_database  = '../data/your_times_here'                               ! Integer ID of weeks past 1 Jan 2025 to consider
    character(*), parameter                 :: output_file    = '../OUTPUT.data'                                        ! Output file

    integer, parameter                      :: output_unit = 97                                                         ! Unit number for storing the final list of candidates

    character(*), parameter                 :: planetary_ephemeris      = '../data/de414.bsp'
    character(*), parameter                 :: leapseconds_file         = '../data/naif0008.tls'
    character(*), parameter                 :: candidate_ephem_prefix   = '../../ephemeris_files/'                      ! Prefix for directory path to candidate ephemerides

    double precision, pointer               :: manifold_data_elements(:,:)                                              ! num_targets x 6

    double precision                        :: manifold_data_state(6)                                                   ! [a, e, i, O, o, M, t] of the target states
    double precision                        :: asteroid_elements(8)                                                     ! Orbital elements of the asteroid
    double precision                        :: asteroid_state(6)                                                        ! State of the asteroid

    integer, allocatable                    :: dataset(:)                                                               ! Holds the weeks to compute
    integer, allocatable                    :: asteroid_database(:)                                                     ! Holds the candidates to compute

end module problem_parameters
