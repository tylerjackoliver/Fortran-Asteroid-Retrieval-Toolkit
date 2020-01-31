PROGRAM MAIN

    ! ///////////////////////////////////////////////////////////
    !
    ! Interfaces the GLOBAL global optimsation routine with the
    ! Lambert transfer cost function
    !
    ! ///////////////////////////////////////////////////////////

    use global_minimum                                                              ! Use the optimiser
    use state_determination                                                         ! Use the state determination routines
    use problem_parameters                                                          ! Problem set-up parameters

    integer :: i
    real*8  :: min_vel

    INTERFACE
        SUBROUTINE funct(x, f, nparm, m)
            IMPLICIT NONE
            DOUBLE PRECISION, INTENT(IN)   :: x(:)
            DOUBLE PRECISION, INTENT(OUT)  :: f
            INTEGER, INTENT(IN)            :: nparm, m
        END SUBROUTINE funct
    END INTERFACE

    ! Load variables into memory

    call MPI_VARIABLE_INIT()

    MIN(1) = 0.90614223E+09
    MIN(2) = 0.10611839E+09

    ! Call the cost function 5 times

    do i = 1,5

        call FUNCT(MIN, MIN_VEL, 2, 1)

    end do 

    ! Exit gracefully

    call VARIABLE_DESTRUCT()

END

    SUBROUTINE FUNCT(X, MIN_VEL, NPARMTR, MM)

        ! ///////////////////////////////////////////////////////////
        !
        ! Optimiser cost function
        !
        ! Inputs
        ! ~~~~~~
        ! X: State vector. X = (epoch, transfer time)
        ! NPARM: Optimiser variable.
        ! M: Optimiser variable.
        !
        ! InOuts
        ! ~~~~~~
        ! MIN_VEL: minimum-velocity transfer in the dataset
        !
        ! ///////////////////////////////////////////////////////////

        use precision_kinds                                                         ! Defines appropriate variable kinds for double and quadruple precision
        use constants                                                               ! Standardises constants
        use state_determination                                                     ! Use the state determination routines (self-written)
        use fortran_astrodynamics_toolkit                                           ! Use the FAT, third-party astrodynamics library
        use problem_parameters                                                      ! Problem set-up parameters

        implicit none

        integer, intent(in)                         :: NPARMTR, MM
            
        real*8, intent(in)                          :: x(:)                         ! Input state vector
        real*8                                      :: state_can(6)                 ! Asteroid candidate state
        real*8		                                :: state_targ(6)                ! Un-rotated target state
        real*8		                                :: state_rot(6)                 ! Rotated target state
        real*8		                                :: state_rot_mod(6)             ! Rotated target state
        real*8                                      :: transfer_epoch               ! Epoch of transfer
        real*8                                      :: transfer_time                ! Time-of-flight
        real*8		                                :: tt 							! Transfer time
        real*8		                                :: vx1, vx2, vy1, vy2, vz1, vz2 ! Velocities of the beginning and end of the Lambert arc
        real*8		                                :: vtx, vty, vtz, vcx, vcy, vcz ! Velocities of the candidate and target at beginning and end of ""
        real*8	                                    :: transfer_v1, transfer_v2     ! Departure velocity of Lambert transfer, insertion velocity of Lambert transfer
        real*8		                                :: transfer_vel                 ! Transfer velocities
        real*8, intent(out)                         :: min_vel                      ! Minimum velocity transfer
        real*8                                      :: v1RET(3)
        real*8                                      :: v2RET(3)
        real*8                                      :: xOut(3)

        logical                                     :: long_way						! Which direction for the Lambert transfer
        logical                                     :: run_ok                       ! Boolean success variable for the Lambert

        real*8, allocatable, dimension(:, :)        :: v1                           ! Velocity of beginning of Lambert arc
        real*8, allocatable, dimension(:, :)        :: v2                           ! Velocity of the end of the Lambert arc

        integer		                                :: multi_rev = 4                ! Number of Lambert arc revolutions (up to)
        integer 	                                :: num_rows                     ! Size of the v1, v2 arrays
        integer 	                                :: i, j, k                      ! Loop sentinels
        integer                                     :: iostate                      ! Status flag for Fortran file I/O
        integer                                     :: best_index                   ! Best index of target state
        integer                                     :: itercount = 0                ! Iteration counter

        ! Initialise the transfer epoch

        transfer_epoch = x(1)

        ! Initialise the transfer time

        tt = x(2)

        ! Initialise velocities

        min_vel = 1.d6

        ! Compute the candidate position

        call CANDIDATE_POSITION(transfer_epoch,state_epoch,state_can_orig,state_can)

        ! Main iteration loop: go through the data file and compute Lamberts to that state

        main_loop: do itercount = 1, num_targets, 2

            ! Get states

            state_targ = dataset(itercount, :)

            ! Rotate the state above via the relations in sanchez et. al.

            call ROTATOR_MOD(state_targ, transfer_epoch, tt, state_rot)

            ! Define initial and final velocities of the target and of
            ! the candidate

            vcx = state_can(4)                                                                                                  ! Candidate; from arguments
            vcy = state_can(5)                                                                                                  ! Candidate; from arguments
            vcz = state_can(6)                                                                                                  ! Candidate; from arguments

            vtx = state_rot(4)                                                                                                  ! Rotated state: from file
            vty = state_rot(5)                                                                                                  ! Rotated state: from file
            vtz = state_rot(6)                                                                                                  ! Rotated state: from file

            ! Now, call the lambert solver on state_rot (the target
            ! state rotated to the correct epoch). First, with long_way
            ! set to false.
            !
            ! TODO: Make this a function

            long_way = .false.

            call solve_lambert_izzo(state_can(1:3), state_rot(1:3), tt, mu, long_way, multi_rev, v1, v2, run_ok)

            if (run_ok .eqv. .false.) then

                write(*,*) "Warning: Izzo's solver failed."

                ! We will re-run the Lambert problem using Gooding's
                ! algorithm, which is typically more robust (at the
                ! expense of speed.)

                call solve_lambert_gooding(state_can(1:3),state_rot(1:3),tt,mu,long_way,multi_rev,v1,v2,run_ok)

                if (run_ok .eqv. .false.) then

                    min_vel = 1e6

                end if

            end if

            ! Determine the number of solutions; we know the array is n x 3, so get the size here

            num_rows = size(v1, 2) ! // Assumes always bigger than 3 //

            do j = 1, num_rows

                vx1 = v1(1,j)                                                                                                   ! x-velocity at start of Lambert arc
                vy1 = v1(2,j)                                                                                                   ! y-velocity at start of Lambert arc
                vz1 = v1(3,j)                                                                                                   ! z-velocity at start of Lambert arc

                vx2 = v2(1,j)                                                                                                   ! x-velocity at end of Lambert arc
                vy2 = v2(2,j)                                                                                                   ! y-velocity at end of Lambert arc
                vz2 = v2(3,j)                                                                                                   ! z-velocity at end of Lambert arc

                transfer_v1  = ((vx1-vcx)**2.d0+(vy1-vcy)**2.d0+(vz1-vcz)**2.d0)**.5d0                                          ! deltaV from the candidate and the beginning of the lambert arc
                transfer_v2  = ((vx2-vtx)**2.d0+(vy2-vty)**2.d0+(vz2-vtz)**2.d0)**.5d0                                          ! deltaV from the target and the end of the lambert arc
                transfer_vel = transfer_v1 + transfer_v2                                                                        ! Total delta v is the sum of both; both > 0

                if (transfer_vel < min_vel) then

                    min_vel = transfer_vel
                    best_index = itercount
                    v1RET = v1(:,j)
                    v2RET = v2(:,j)
                    xOut = state_rot(1:3)

                end if

            end do

            ! Try all of the above, but with the long-way solution
            ! This should capture cases where the dV signs cancel, in case
            ! 

            long_way = .true.

            call solve_lambert_izzo(state_can(1:3),state_rot(1:3),tt, mu, long_way,multi_rev,v1,v2,run_ok)

            if (run_ok .eqv. .false.) then

                write(*,*) "Warning: Izzo's solver failed."

                ! We will re-run the Lambert problem using Gooding's
                ! algorithm, which is typically more robust (at the
                ! expense of speed.)

                call solve_lambert_gooding(state_can(1:3),state_rot(1:3),tt,mu,long_way,multi_rev,v1,v2,&
                    run_ok)
                    
                if (run_ok .eqv. .false.) then

                    min_vel = 1e6

                end if

            end if

            ! Determine the number of solutions; we know the array is n x 3, so get the size here

            num_rows = size(v1, 2)

            do j = 1, num_rows

                vx1 = v1(1,j)                                                                                                   ! x-velocity at start of Lambert arc
                vy1 = v1(2,j)                                                                                                   ! y-velocity at start of Lambert arc
                vz1 = v1(3,j)                                                                                                   ! z-velocity at start of Lambert arc

                vx2 = v2(1,j)                                                                                                   ! x-velocity at end of Lambert arc
                vy2 = v2(2,j)                                                                                                   ! y-velocity at end of Lambert arc
                vz2 = v2(3,j)                                                                                                   ! z-velocity at end of Lambert arc

                transfer_v1  = ((vx1-vcx)**2.d0+(vy1-vcy)**2.d0+(vz1-vcz)**2.d0)**.5d0                                          ! deltaV from the candidate and the beginning of the lambert arc
                transfer_v2  = ((vx2-vtx)**2.d0+(vy2-vty)**2.d0+(vz2-vtz)**2.d0)**.5d0                                          ! deltaV from the target and the end of the lambert arc
                transfer_vel = transfer_v1 + transfer_v2                                                                        ! Total delta v is the sum of both; both > 0

                if (transfer_vel < min_vel) then

                    min_vel  = transfer_vel
                    best_index = itercount
                    v1RET = v1(:,j)
                    v2RET = v2(:,j)
                    xOut = state_rot(1:3)

                end if

            end do

    end do main_loop

    rewind(69)
    close(69)

    print *, "Function evaluation completed. Top transfer was"
    print *, "Min vel: ", min_vel
    print *, "Best index: ", best_index
    print *, "v1RET: ", v1RET
    print *, "v2RET: ", v2RET
    print *, "xOut: ", xOut

    RETURN

END
