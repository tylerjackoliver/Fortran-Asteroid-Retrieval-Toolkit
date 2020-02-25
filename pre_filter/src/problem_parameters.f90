module problem_parameters

    use precision
    use constants

    implicit none

    character(*), parameter                 :: target_dataset = '../data/piBy8Targets.csv'
    character(*), parameter                 :: spk_database   = '../data/DatabaseSPKs.csv'
    character(*), parameter                 :: output_file    = '../OUTPUT.data'

    integer, parameter                      :: output_unit = 97                                           ! Unit number for storing the final list of candidates

    character(*), parameter                 :: planetary_ephemeris      = '../data/de414.bsp'
    character(*), parameter                 :: leapseconds_file         = '../data/naif0008.tls'
    character(*), parameter                 :: candidate_ephem_prefix   = '../data/ephem/'

    real(kind=dp), allocatable              :: target_elements(:,:)                                  ! [a, e, i, O, o, M, t] of the target states
    real(kind=dp), allocatable              :: manifold_data_elements(:,:)                           ! num_targets x 6
    real(kind=dp), allocatable              :: asteroid_elements(8)                                  ! Orbital elements of the asteroid
    real(kind=dp), allocatable              :: asteroid_state(6)                                     ! State of the asteroid

end module problem_parameters