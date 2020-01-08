module lambert_transfer_module

    use precision_kinds
    use constants
    use fortran_astrodynamics_toolkit
    
    implicit none

    contains

        subroutine lambert_transfer(state_can, state_rot, tt, multi_rev, min_vel)

            real*8,  intent(in)                         :: state_can(6)                 ! Asteroid candidate state
            real*8,  intent(in)                         :: state_rot(6)                 ! Rotated target state
            real*8,  intent(in)                         :: tt 							! Transfer time

            integer, intent(in)                         :: multi_rev                    ! Number of Lambert revolutions to consider

            real*8,  intent(out)                        :: min_vel                      ! Minimum velocity transfer found

            real*8		                                :: state_targ(6)                ! Un-rotated target state
            real*8                                      :: transfer_epoch               ! Epoch of transfer
            real*8                                      :: transfer_time                ! Time-of-flight
            real*8		                                :: vx1, vx2, vy1, vy2, vz1, vz2 ! Velocities of the beginning and end of the Lambert arc
            real*8		                                :: vtx, vty, vtz, vcx, vcy, vcz ! Velocities of the candidate and target at beginning and end of ""
            real*8	                                    :: transfer_v1, transfer_v2     ! Departure velocity of Lambert transfer, insertion velocity of Lambert transfer
            real*8		                                :: transfer_vel                 ! Transfer velocities

            real*8, allocatable, dimension(:, :)        :: v1                           ! Velocity of beginning of Lambert arc
            real*8, allocatable, dimension(:, :)        :: v2                           ! Velocity of the end of the Lambert arc

            integer                                     :: num_rows                     ! Number of rows in the Lambert transfer array
            integer                                     :: i, j, k

            logical                                     :: long_way
            logical                                     :: run_ok                       ! Lambert transfer success flag

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

            long_way = .false.

            ! Initialise minimum velocity

            min_vel = 1.e6

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

                end if

            end do

        end subroutine lambert_transfer

end module lambert_transfer_module