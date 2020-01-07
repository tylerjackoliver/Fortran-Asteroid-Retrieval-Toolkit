PROGRAM MAIN

    ! ///////////////////////////////////////////////////////////
    !
    ! Interfaces the GLOBAL global optimsation routine with the
    ! Lambert transfer cost function
    !
    ! ///////////////////////////////////////////////////////////

    USE global_minimum                                                              ! Use the optimiser
    use state_determination                                                         ! Use the state determination routines

    real*8  X0(30,20), F0(20), MIN(2), MAX(2), TRANSFER_EPOCH, VALUE                ! Optimiser variables
    real*8  state_can_orig(6)                                                       ! Un-rotated asteroid candidate state
    real*8  time_lower                                                              ! Lower bound for the optimiser (epoch)
    real*8  time_upper                                                              ! Upper bound for the optimiser (epoch)
    real*8  state_epoch                                                             ! Epoch of the candidate state
    
    integer M, NPARM                                                                ! More optimiser variables

    M      = 1
    NPARM  = 2
    NSAMPL = 50
    NSEL   = 2
    IPR    = 77

    ! Open the output file for the optimiser

    OPEN(IPR, FILE='OUTPUT')

    ! Number of significant figures in optimiser print-outs

    NSIG = 9

    ! Initialise transfer epoch bounds - call state_finder to get bounds and ignore
    ! other outputs (state_can_orig, state_epoch)

    call STATE_FINDER('3550232', state_can_orig, state_epoch, time_lower, time_upper)

    ! Initialise bounds for the optimiser

    MIN(1) = time_lower                                                             ! Minimum transfer epoch (ephemeris seconds)
    MIN(2) = 0.0D0 * 86400.D0                                                       ! Minimum transfer duration (seconds)
    MAX(1) = time_upper                                                             ! Maximum transfer epoch (ephemeris seconds)
    MAX(2) = 1600.D0 * 86400.D0                                                     ! Maximum transfer epoch (seconds)

    ! Call the optimiser

    CALL GLOBAL(MIN, MAX, NPARM, M, NSAMPL, NSEL, IPR, NSIG, X0, NC, F0)

    ! Close the output file

    CLOSE(IPR)

END

    SUBROUTINE FUNCT(X, MIN_VEL, NPARM, M)

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

        implicit none

        integer		                                :: NPARM, M
            
        real*8		                                :: x(:)                         ! Input state vector
        real*8                                      :: state_can(6)                 ! Asteroid candidate state
        real*8                                      :: state_can_orig(6)            ! Un-rotated asteroid candidate state
        real*8		                                :: state_targ(6)                ! Un-rotated target state
        real*8		                                :: state_rot(6)                 ! Rotated target state
        real*8                                      :: transfer_epoch               ! Epoch of transfer
        real*8                                      :: transfer_time                ! Time-of-flight
        real*8                                      :: state_epoch                  ! Epoch of the candidate state
        real*8		                                :: time_lower                   ! Lower bound for the optimisation time - dummy here
        real*8		                                :: time_upper                   ! Upper bound for the optimisation time - dummy here
        real*8		                                :: tt 							! Transfer time
        real*8		                                :: vx1, vx2, vy1, vy2, vz1, vz2 ! Velocities of the beginning and end of the Lambert arc
        real*8		                                :: vtx, vty, vtz, vcx, vcy, vcz ! Velocities of the candidate and target at beginning and end of ""
        real*8	                                    :: transfer_v1, transfer_v2     ! Departure velocity of Lambert transfer, insertion velocity of Lambert transfer
        real*8		                                :: transfer_vel                 ! Transfer velocities
        real*8		                                :: min_vel                      ! Minimum velocity transfer

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

        character(*), parameter                     :: targ_can='3550232'           ! Target candidate string
        
        ! Find the state of the target and the original
    
        call STATE_FINDER(targ_can, state_can_orig, state_epoch, time_lower, time_upper)

        ! Get the candidate position at the correct time x(1)

        transfer_epoch = x(1)

        ! Initialise the transfer time

        tt = x(2)

        ! Initialise velocities

        min_vel = 1.d6

        ! Open input file

        open(69, file='../data/2019-11-20_L2PlanarBackCondsGlobal.csv')

        ! Compute the candidate position

        call CANDIDATE_POSITION(transfer_epoch,state_epoch,state_can_orig,state_can)

        itercount = 0

        ! Main iteration loop: go through the data file and compute Lamberts to that state

        main_loop: do

            ! Read state in

            read(69, *, iostat=iostate) state_targ(1), state_targ(2), state_targ(3), state_targ(4), state_targ(5), state_targ(6)

            ! If at the end of the file, exit gracefully

            if (iostate .ne. 0) then

                exit main_loop
    
            end if

            ! In case the optimiser shits the bed with invalid inputs

            if (tt .lt. 0) then

                exit main_loop

            end if

            ! Increase the iteration counter by 1

            itercount = itercount + 1

            ! Rotate the state above via the relations in sanchez et. al.

            call ROTATOR(state_targ, transfer_epoch, tt, state_rot)

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

                end if

        end do

    end do main_loop

    rewind(69)
    close(69)

    RETURN

    END
