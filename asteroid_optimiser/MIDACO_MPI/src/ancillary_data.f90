module ancillary_data

    use precision_kinds
    use constants
    use problem_parameters
    use variable_initialisation
    use fortran_astrodynamics_toolkit
    use compute_spline

    implicit none

    contains

        subroutine generate_ancillary_data(t_epoch, transfer_time, t_end, n_mnfd, orbit_num, &
            r1, v1O, r2, v2O, tmani, target_point)

            implicit none

            double precision, intent(in)                :: t_epoch
            double precision, intent(in)                :: transfer_time
            double precision, intent(in)                :: t_end
            double precision, intent(in)                :: n_mnfd
            
            integer,          intent(out)               :: orbit_num

            double precision, intent(out)               :: r1(6)
            double precision, intent(out)               :: v1O(3)
            double precision, intent(out)               :: r2(6)
            double precision, intent(out)               :: v2O(3)
            double precision, intent(out)               :: target_point(6)
            double precision, intent(out)               :: tmani

            integer                                     :: i, j

            double precision                            :: state_can(6)                 ! Asteroid candidate state
            double precision		                    :: state_targ(6)                ! Un-rotated target state
            double precision		                    :: state_rot(6)                 ! Rotated target state
            double precision		                    :: vx1, vx2, vy1, vy2, vz1, vz2 ! Velocities of the beginning and end of the Lambert arc
            double precision		                    :: vtx, vty, vtz, vcx, vcy, vcz ! Velocities of the candidate and target at beginning and end of ""
            double precision	                        :: transfer_v1, transfer_v2     ! Departure velocity of Lambert transfer, insertion velocity of Lambert transfer
            double precision		                    :: transfer_vel                 ! Transfer velocities
            double precision		                    :: min_vel                      ! Minimum velocity transfer
    
            logical                                     :: long_way						! Which direction for the Lambert transfer
            logical                                     :: run_ok                       ! Boolean success variable for the Lambert
    
            double precision, allocatable, dimension(:, :)        :: v1                           ! Velocity of beginning of Lambert arc
            double precision, allocatable, dimension(:, :)        :: v2                           ! Velocity of the end of the Lambert arc
    
            integer		                                :: multi_rev = 4                ! Number of Lambert arc revolutions (up to)
            integer 	                                :: num_rows                     ! Size of the v1, v2 arrays
            integer 	                                :: j                            ! Loop sentinels
            integer                                     :: best_index                   ! Best index of target state
            integer                                     :: itercount = 0                ! Iteration counter
            integer                                     :: orbit_choice                 ! The orbit we're using - bruting over this
    
            ! Initialise velocities
    
            min_vel = 1.d6
    
            ! Compute the candidate position
    
            call CANDIDATE_POSITION(t_epoch, state_can)

            r1 = state_can
    
            ! Main iteration loop: go through the data file and compute Lamberts to that state
    
            main_loop: do orbit_choice = 1, num_orbits ! Brute over J
    
                ! Get states via interpolation of the dataset
    
                call BSPLINE_INTERPOLATE(t_end, n_mnfd, orbit_choice, state_targ)
    
                ! Rotate the state above via the relations in sanchez et. al.
    
                call ROTATOR(state_targ, t_epoch, transfer_time, state_rot)
    
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
    
                call solve_lambert_izzo(state_can(1:3), state_rot(1:3), transfer_time, mu, long_way, multi_rev, v1, v2, run_ok)
    
                if (run_ok .eqv. .false.) then
    
                    write(*,*) "Warning: Izzo's solver failed."
    
                    ! We will re-run the Lambert problem using Gooding's
                    ! algorithm, which is typically more robust (at the
                    ! expense of speed.)
    
                    call solve_lambert_gooding(state_can(1:3),state_rot(1:3),transfer_time,mu,long_way,multi_rev,v1,v2,run_ok)
    
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
                        orbit_num = orbit_choice
                        r2 = state_rot
                        v1O = (/(vx1-vcx), (vy1-vcy), (vz1-vcz)/)
                        v2O = (/(vx2-vtx), (vy2-vty), (vz2-vtz)/)
                        
                    end if
    
                end do
    
                ! Try all of the above, but with the long-way solution
                ! This should capture cases where the dV signs cancel, in case
    
                long_way = .true.
    
                call solve_lambert_izzo(state_can(1:3),state_rot(1:3),transfer_time, mu, long_way,multi_rev,v1,v2,run_ok)
    
                if (run_ok .eqv. .false.) then
    
                    write(*,*) "Warning: Izzo's solver failed."
    
                    ! We will re-run the Lambert problem using Gooding's
                    ! algorithm, which is typically more robust (at the
                    ! expense of speed.)
    
                    call solve_lambert_gooding(state_can(1:3),state_rot(1:3),transfer_time,mu,long_way,multi_rev,v1,v2,&
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
    
                        min_vel = transfer_vel
                        orbit_num = orbit_choice
                        r2 = state_rot
                        v1O = (/(vx1-vcx), (vy1-vcy), (vz1-vcz)/)
                        v2O = (/(vx2-vtx), (vy2-vty), (vz2-vtz)/)
                        
                    end if
    
                end do
    
        end do main_loop

        call generate_target_point(n_mnfd, J, target_point, tmani)

    end subroutine generate_ancillary_data

    subroutine generate_target_point(n_mnfd, J, target_point, tmani)

        use precision
        use constants
        use problem_parameters
        use variable_initialisation
        use compute_spline

        implicit none

        real(kind=dp), intent(in)               :: n_mnfd
        integer,       intent(in)               :: J

        real(kind=dp), intent(out)              :: target_point(6)
        real(kind=dp), intent(out)              :: tmani

        call b_spline_interpolate_perturbed_conds(n_mnfd, J, target_point, tmani)

    end subroutine generate_target_point

end module ancillary_data