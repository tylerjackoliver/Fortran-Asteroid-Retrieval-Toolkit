PROGRAM MAIN

    ! ///////////////////////////////////////////////////////////
    !
    ! Interfaces the GLOBAL global optimsation routine with the
    ! Lambert transfer cost function
    !
    ! ///////////////////////////////////////////////////////////

    use state_determination                                                         ! Use the state determination routines
    use problem_parameters                                                          ! Problem set-up parameters
    use variable_initialisation
    use ancillary_data

    implicit none

    integer      :: target_count
    character(7) :: targ_can_temp
    
    logical      :: first_load = .true.

    ! Initialise MPI functionality; read-in datasets

    call MPI_VARIABLE_INIT()
    
    ! Perform the chunk of our work for the dataset

    do target_count = (mpi_id_world+1), size(targ_can_array), mpi_world_size

        write(targ_can_temp, '(I7)'), targ_can_array(target_count)
        targ_can = trim(targ_can_temp)

        if (first_load) then

            call variable_init()
            first_load = .false.    

        else

            call intermediate_variable_init()

        end if

        call run_optim()
        call get_pareto_front()
        call intermediate_variable_destruct()

    end do

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    ! Exit gracefully

    ! call MPI_VARIABLE_DESTRUCT()

contains

    subroutine run_optim()

        ! Print header

        call midaco_print(1, print_eval, save_to_file, optim_flag, optim_stop, F, G, XOPT, &
                        XL, XU, O, N, NI, M, ME, RW, PF, max_eval, max_time, param, 1, 0, key)

        do while (optim_stop .eq. 0) 

            ! Evaluate objective function and constraints (none)

            call problem_function(f, xopt)

            ! Call MIDACO

            call midaco(1, O, N, NI, M, ME, XOPT, F, G, XL, XU, optim_flag, optim_stop, &
                        param, rw, lrw, iw, liw, pf, lpf, key)

            ! Print again

            call midaco_print(2, print_eval, save_to_file, optim_flag, optim_stop, F, G, XOPT, &
            XL, XU, O, N, NI, M, ME, RW, PF, max_eval, max_time, param, 1, 0, key)

    end do

    end subroutine run_optim

    subroutine problem_function(F, X)

        ! Output: F, given contraints G and X

        use precision_kinds

        implicit none

        double precision, intent(out)   :: F(2)
        double precision, intent(in)    :: X(4)

        double precision :: min_vel

        call funct(x, min_vel) ! deltaV
        F(1) = min_vel
        F(2) = x(2)            ! tt - ignored, only for Pareto front

    end subroutine problem_function

    SUBROUTINE FUNCT(X, MIN_VEL)

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
        use compute_spline
        use variable_initialisation

        implicit none

        double precision		                                :: x(4)                         ! Input state vector
        double precision                                        :: state_can(6)                 ! Asteroid candidate state
        double precision		                                :: state_targ(6)                ! Un-rotated target state
        double precision                                        :: state_dim(6)
        double precision		                                :: state_rot(6)                 ! Rotated target state
        double precision                                        :: transfer_epoch               ! Epoch of transfer
        double precision		                                :: tt 							! Transfer time
        double precision                                        :: t_end                        ! Desired backwards integration time
        double precision                                        :: n_mnfd                       ! Desired point along the orbit
        double precision		                                :: vx1, vx2, vy1, vy2, vz1, vz2 ! Velocities of the beginning and end of the Lambert arc
        double precision		                                :: vtx, vty, vtz, vcx, vcy, vcz ! Velocities of the candidate and target at beginning and end of ""
        double precision	                                    :: transfer_v1, transfer_v2     ! Departure velocity of Lambert transfer, insertion velocity of Lambert transfer
        double precision		                                :: transfer_vel                 ! Transfer velocities
        double precision		                                :: min_vel                      ! Minimum velocity transfer

        logical                                                 :: long_way						! Which direction for the Lambert transfer
        logical                                                 :: run_ok                       ! Boolean success variable for the Lambert

        double precision, allocatable, dimension(:, :)          :: v1                           ! Velocity of beginning of Lambert arc
        double precision, allocatable, dimension(:, :)          :: v2                           ! Velocity of the end of the Lambert arc

        integer		                                            :: multi_rev = 4                ! Number of Lambert arc revolutions (up to)
        integer 	                                            :: num_rows                     ! Size of the v1, v2 arrays
        integer 	                                            :: j                            ! Loop sentinels
        integer                                                 :: best_index                   ! Best index of target state
        integer                                                 :: itercount = 0                ! Iteration counter
        integer                                                 :: orbit_choice                 ! The orbit we're using - bruting over this

        ! Initialise variables from input vector

        transfer_epoch  = x(1)
        tt              = x(2)
        t_end           = x(3)
        n_mnfd          = x(4)

        if (x(1) .lt. 0 .or. x(2) .lt. 0 .or. x(3) .lt. 0 .or. x(4) .lt. 0) then

            min_vel = 1.d6
            RETURN

        end if

        ! Initialise velocities

        min_vel = 1.d6

        ! Compute the candidate position

        call CANDIDATE_POSITION(transfer_epoch, state_can)

        ! Main iteration loop: go through the data file and compute Lamberts to that state

        main_loop: do orbit_choice = 1, num_orbits ! Brute over J

            ! Get states via interpolation of the dataset

            call BSPLINE_INTERPOLATE(t_end, n_mnfd, orbit_choice, state_targ)

            ! Rotate into the synodic frame

            call GLOBAL_ROTATE(state_targ, transfer_epoch, state_dim)

            ! Dimensionalise

            state_dim(1:3) = state_dim(1:3) * au
            state_dim(4:6) = state_dim(4:6) * au * 2.d0 * 4.d0 * datan(1.d0)/ (86400.d0 * 365.25d0)

            ! Rotate the state above via the relations in sanchez et. al.

            call ROTATOR(state_dim, transfer_epoch, tt, state_rot)

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

            ! print *, state_rot

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

    end do main_loopmem_free

    RETURN

END

end program main
