program main

    use constants
    use lambert_module

    implicit none

    character(*), parameter :: input_file='../data/planarL2_newRange_globalBackwards.dat'

    double precision :: state_targ(6)

    double precision :: state_can(6)

    double precision :: tt = 1239.d0 * 86400

    REAL(kind=dp)                              :: vx1                 ! x-component of velocity at beginning of Lambert arc
    REAL(kind=dp)                              :: vy1                 ! y-component of velocity at beginning of Lambert arc
    REAL(kind=dp)                              :: vz1                 ! z-component of velocity at beginning of Lambert arc
    REAL(kind=dp)                              :: vx2                 ! x-component of velocity at end of Lambert arc
    REAL(kind=dp)                              :: vy2                 ! y-component of velocity at end of Lambert arc
    REAL(kind=dp)                              :: vz2                 ! z-component of velocity at end of Lambert arc
    REAL(kind=dp)                              :: transfer_v1         ! Transfer velocity required at start of arc
    REAL(kind=dp)                              :: transfer_v2         ! Transfer velocity required at end of arc
    REAL(kind=dp)                              :: transfer_vel        ! Total transfer velocity
    REAL(kind=dp)                              :: vcx                 ! x-velocity of the candidate
    REAL(kind=dp)                              :: vcy                 ! y-velocity of the candidate
    REAL(kind=dp)                              :: vcz                 ! z-velocity of the candidate
    REAL(kind=dp)                              :: vtx                 ! x-velocity of the target
    REAL(kind=dp)                              :: vty                 ! y-velocity of the target
    REAL(kind=dp)                              :: vtz                 ! z-velocity of the target

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: v1                  ! ALLOCATABLE attribute required by the Fortran Astrodynamics Toolkit
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: v2                  ! ALLOCATABLE attribute required by the Fortran Astrodynamics Toolkit

    INTEGER(kind=sp)                           :: i, j, num_rows, multi_rev
    logical                                    :: long_way, run_ok

    open(unit=80, file=input_file)

    read(80,*) state_targ(1), state_targ(2), state_targ(3), state_targ(4), state_targ(5), state_targ(6)

                                         

    state_can(1) = 19609904.356833760
    state_can(2) = 155271474.52677387
    state_can(3) = -480329.37284089625

    state_can(4) = -28.669817588147247
    state_can(5) = 5.3403701099956571
    state_can(6) = -8.3903904761525316E-002

    vcx = state_can(4)                                                                                                  ! Candidate; from arguments
    vcy = state_can(5)                                                                                                  ! Candidate; from arguments
    vcz = state_can(6)                                                                                                  ! Candidate; from arguments

    vtx = state_targ(4)                                                                                                  ! Rotated state: from file
    vty = state_targ(5)                                                                                                  ! Rotated state: from file
    vtz = state_targ(6)                                                                                                  ! Rotated state: from file

    ! Now, call the lambert solver on state_rot (the target
    ! state rotated to the correct epoch). First, with long_way
    ! set to false

    multi_rev = 3

    do i = 1, 100000

        long_way = .false.

        call solve_lambert_izzo(state_can(1:3), state_targ(1:3), tt, mu, long_way, multi_rev, v1, v2, run_ok)

        ! Determine the number of solutions; we know the array is n x 3, so get the size here

        num_rows = size(v1) ! // Assumes always bigger than 3 //

        do j = 1, num_rows

            vx1 = v1(1,j)                                                                                                   ! x-velocity at start of Lambert arc
            vy1 = v1(2,j)                                                                                                   ! y-velocity at start of Lambert arc
            vz1 = v1(3,j)                                                                                                   ! z-velocity at start of Lambert arc

            vx2 = v2(1,j)                                                                                                   ! x-velocity at end of Lambert arc
            vy2 = v2(2,j)                                                                                                   ! y-velocity at end of Lambert arc
            vz2 = v2(3,j)                                                                                                   ! z-velocity at end of Lambert arc

            transfer_v1  = ((vx1-vcx)**2+(vy1-vcy)**2+(vz1-vcz)**2)**.5                                                     ! deltaV from the candidate and the beginning of the lambert arc
            transfer_v2  = ((vx2-vtx)**2+(vy2-vty)**2+(vz2-vtz)**2)**.5                                                     ! deltaV from the target and the end of the lambert arc
            transfer_vel = transfer_v1 + transfer_v2                                                                        ! Total delta v is the sum of both; both > 0

        end do

        write(*,*) "Transfer vel is", transfer_vel

    end do

    close(80)

end program main
