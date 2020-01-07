SUBROUTINE CALFUN (N,X,F)
        
        use precision_kinds                                 ! Defines appropriate variable kinds for double and quadruple precision
        use constants                                       ! Standardises constants
        use state_determination                             ! Use the state determination routines (self-written)
        use fortran_astrodynamics_toolkit                   ! Use the FAT, third-party astrodynamics library

        implicit none

        integer		                        :: N
        
        real*8		                        :: x(n), f                              ! Input state vector (x = (epoch, transfer_time)
        real*8                              :: state_can(6)                         ! Asteroid candidate state
        real*8                              :: state_can_orig(6)                    ! Un-rotated asteroid candidate state
        real*8		                        :: state_targ(6)                        ! State of the target, un-rotated
        real*8		                        :: state_rot(6)                         ! State of the target, rotated
        real*8                              :: transfer_epoch                       ! Epoch of transfer
        real*8                              :: transfer_time                        ! Time-of-flight
        real*8                              :: state_epoch                          ! Epoch of the candidate state
        real*8		                        :: time_lower                           ! Lower bound of the 
        real*8		                        :: time_upper
        real*8		                        :: tt 									! Transfer time
        real*8		                        :: vx1, vx2, vy1, vy2, vz1, vz2
        real*8		                        :: vtx, vty, vtz, vcx, vcy, vcz
        real*8	                            :: transfer_v1, transfer_v2
        real*8		                        :: transfer_vel
        real*8		                        :: min_vel = 10.d6

        logical                             :: long_way								
        logical                             :: run_ok

        real*8, allocatable, dimension(:, :):: v1
        real*8, allocatable, dimension(:, :):: v2

        integer		                        :: multi_rev = 4
        integer 	                        :: num_rows
        integer 	                        :: i, j
        integer                             :: iostate

        character(*), parameter  :: targ_can='3550232'                              ! Target candidate string
        
        ! Find the state of the target and the original
    
        call STATE_FINDER(targ_can, state_can_orig, state_epoch, time_lower, time_upper)

        ! Get the candidate position at the correct time x(1)

        transfer_epoch = x(1)

        ! Initialise the transfer time

        tt = x(2)

        ! Open input file

        open(69, file='../data/2019-11-20_L2PlanarBackCondsGlobal.csv')

        ! Compute the candidate position

        call CANDIDATE_POSITION(transfer_epoch,state_epoch,state_can_orig,state_can)

        ! ! Initialise the target state at J2000

        ! state_targ = (/-141002146.07237700, 87602033.928656995, 0.00, -13.8609380000, -24.0032210000, 0.00/)

        do

            read(69, *, iostat=iostate) state_targ(1), state_targ(2), state_targ(3), state_targ(4), state_targ(5), state_targ(6)

            if (iostate .ne. 0) then

                exit
    
            end if

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
            ! set to false

            long_way = .false.

            call solve_lambert_izzo(state_can(1:3), state_rot(1:3), tt, mu, long_way, multi_rev, v1, v2, run_ok)

            if (run_ok .eqv. .false.) then

                write(*,*) "Warning: Izzo's solver failed."

                ! We will re-run the Lambert problem using Gooding's
                ! algorithm, which is typically more robust (at the
                ! expense of speed.)

                call solve_lambert_gooding(state_can(1:3),state_rot(1:3),tt,mu,long_way,multi_rev,v1,v2,run_ok)

                if (run_ok .eqv. .false.) then

                    error stop "Lambert methods failed. Aborting..."

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

                    error stop "Lambert methods failed. Aborting..."

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

                transfer_v1  = ((vx1-vcx)**2.d0+(vy1-vcy)**2.d0+(vz1-vcz)**2.d0)**.5d0                                                     ! deltaV from the candidate and the beginning of the lambert arc
                transfer_v2  = ((vx2-vtx)**2.d0+(vy2-vty)**2.d0+(vz2-vtz)**2.d0)**.5d0                                                      ! deltaV from the target and the end of the lambert arc
                transfer_vel = transfer_v1 + transfer_v2                                                                        ! Total delta v is the sum of both; both > 0

                if (transfer_vel < min_vel) then

                    min_vel  = transfer_vel

                end if

        end do

    end do

    f = min_vel

    close(69)

    RETURN
END
