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
    use midaco_interface

    implicit none

    integer             :: target_time_array
    integer             :: target_count ! For now, fix to be only one candidate at a time
    integer             :: tt_counter
    integer             :: initial_tt

    character(7)        :: targ_can_temp
    character(10)       :: tt_str

    character(4)        :: direction_str
    integer             :: direction, direction_temp
    integer             :: loop_bound    

    logical             :: first_load = .true.

    double precision    :: initial_xopt(4), temp

    ! Initialise MPI functionality; read-in datasets

    call MPI_VARIABLE_INIT()

    !
    ! Open the file with the target candidate solutions in it
    !

    do tt_counter = initial_tt, loop_bound, 3 * direction

        transfer_time = transfer_time_array(tt_counter)

        write(tt_str, '(I4)') tt_counter
        print *, "Time", tt_str

        if (first_load) then

                call variable_init()
                first_load = .false.

        else

                call intermediate_variable_init()

        end if

        ! First run, completely from scratch

        call run_global_optim()

        if (mpi_id_world .eq. 0) then

        call funct(xopt, pareto_front_minimum)
        write(iunit, *) pareto_front_minimum, tt_counter, xopt

        end if
       
        initial_xopt = xopt

        call intermediate_variable_destruct()

    end do
        
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    call MPI_FINALIZE(mpi_err)

    ! call MPI_VARIABLE_DESTRUCT()

contains

    subroutine run_global_optim()

        !! MIDACO-internal parallelism variables

        integer             :: variable_num
        integer             :: thread_num
        integer             :: P                    ! Number of threads
        integer             :: status(MPI_STATUS_SIZE)
        integer             :: MINIMUM_LOCATION

        ! Main results arrays
        
        double precision    :: XOPT_PAR(1000)     ! Increase alloc. if necessary
        double precision    :: F_PAR(100)          ! Increase alloc. if necessary
        double precision    :: G_PAR(100)          ! Increase alloc.

        ! Result/design arrays for slaves

        double precision    :: DX(N), DF(N), DG(N)  ! XOPT, F, G OUT
        double precision    :: EX(N), EF(N), EG(N)  ! XOPT, F, G OUT

        ! For randomly assigning initial conditions

        double precision    :: DIFFERENCE_VECTOR(4) ! FOR XU - XL
        double precision    :: eta                  ! How far along (0 <= x < 1)

        DIFFERENCE_VECTOR = XU - XL
        P = mpi_world_size

        ! If master, initialise MPI and run

        if (mpi_id_world .eq. 0) then

            !
            ! Copy the starting point XOPT into parallel array XOPT_PAR
            !

            do thread_num = 2, P

                do variable_num = 1, N

                    call RANDOM_NUMBER(eta)
                    XOPT_PAR((thread_num-1)*N+variable_num) = XL(variable_num) + &
                        DIFFERENCE_VECTOR(variable_num) * eta

                end do

            end do

            ! Wonderful - now call MPI as master

            do while (optim_stop .eq. 0)

                ! Store variables XOPT in dummy DX, send to slaves

                do thread_num = 2, P

                    do variable_num = 1, N

                        DX(variable_num) = XOPT_PAR((thread_num-1)*N+variable_num)

                    end do

                    ! Blocking send

                    CALL MPI_SEND(DX, N, MPI_DOUBLE_PRECISION, thread_num-1, 1, &
                        MPI_COMM_WORLD, mpi_err)

                end do

                ! Evaluate on master

                call problem_function(F_PAR, XOPT_PAR)

                ! Collect other results

                do thread_num = 2, P

                    call MPI_RECV(DF, O, MPI_DOUBLE_PRECISION, thread_num-1,2,&
                        MPI_COMM_WORLD, status, mpi_err)

                    do variable_num = 1, O

                        F_PAR((thread_num-1)*O+variable_num) = DF(variable_num)

                    end do

                end do

                ! Call MIDACO

                call MIDACO(P, O, N, NI, M, ME, XOPT_PAR, F_PAR, G_PAR, XL, XU, optim_flag, &
                    optim_stop, param, rw, lrw, iw, liw, pf, lpf, save_to_file, max_eval, &
                    max_time, print_eval)
       

                do thread_num = 2, P

                    call MPI_SEND(optim_stop, 1, MPI_INTEGER, thread_num-1,&
                        4, MPI_COMM_WORLD, mpi_err)

                end do

            end do

        else

            optim_stop = 0

            do while (optim_stop .eq. 0)

                call MPI_RECV(EX, N, MPI_DOUBLE_PRECISION, 0, 1, &
                    MPI_COMM_WORLD, STATUS, mpi_err)
                call problem_function(EF, EX)
                
                call MPI_SEND(EF, O, MPI_DOUBLE_PRECISION, 0, 2, &
                    MPI_COMM_WORLD, mpi_err)

                call MPI_RECV(optim_stop, 1, MPI_INTEGER, 0, 4, &
                    MPI_COMM_WORLD, status, mpi_err)

            end do

        end if

        ! Initialise XOPT to best from last

        XOPT = XOPT_PAR(1:N)

    end subroutine run_global_optim

    subroutine problem_function(F, X)

        ! Output: F, given contraints G and X

         

        implicit none

        double precision, intent(out)   :: F(O)
        double precision, intent(in)    :: X(N)

        double precision                :: min_vel

        call funct(x, min_vel) ! deltaV

        F(1) = min_vel

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

                                                                  ! Defines appropriate variable kinds for double and quadruple precision
        use constants                                                               ! Standardises constants
        use state_determination                                                     ! Use the state determination routines (self-written)
        use fortran_astrodynamics_toolkit                                           ! Use the FAT, third-party astrodynamics library
        use problem_parameters                                                      ! Problem set-up parameters
        use compute_spline
        use variable_initialisation

        implicit none

        double precision		                                :: x(N)                         ! Input state vector
        double precision                                        :: state_can(6)                 ! Asteroid candidate state
        double precision		                                :: state_targ(6)                ! Un-rotated target state
        double precision                                        :: state_dim(6)
        double precision		                                :: state_rot(6)                 ! Rotated target state
        double precision                                        :: transfer_epoch               ! Epoch of transfer
        double precision		                                :: tt 							! Transfer time
        double precision                                        :: t_end                        ! Desired backwards integration time
        double precision                                        :: n_mnfd                       ! Desired point along the orbit
        
        double precision                                        :: orbit_choice                 ! The orbit we're using

        double precision		                                :: vx1, vx2, vy1, vy2, vz1, vz2 ! Velocities of the beginning and end of the Lambert arc
        double precision		                                :: vtx, vty, vtz, vcx, vcy, vcz ! Velocities of the candidate and target at beginning and end of ""
        double precision	                                    :: transfer_v1, transfer_v2     ! Departure velocity of Lambert transfer, insertion velocity of Lambert transfer
        double precision		                                :: transfer_vel                 ! Transfer velocities
        double precision		                                :: min_vel                      ! Minimum velocity transfer
        double precision                                        :: upper_limit_time

        logical                                                 :: long_way						! Which direction for the Lambert transfer
        logical                                                 :: run_ok                       ! Boolean success variable for the Lambert

        double precision, allocatable, dimension(:, :)          :: v1                           ! Velocity of beginning of Lambert arc
        double precision, allocatable, dimension(:, :)          :: v2                           ! Velocity of the end of the Lambert arc

        integer		                                            :: multi_rev = 4                ! Number of Lambert arc revolutions (up to)
        integer 	                                            :: num_rows                     ! Size of the v1, v2 arrays
        integer 	                                            :: j                            ! Loop sentinels
        integer                                                 :: best_index                   ! Best index of target state
        integer                                                 :: itercount = 0                ! Iteration counter

        ! Initialise variables from input vector

        transfer_epoch  = x(1)
        t_end           = x(2)
        n_mnfd          = x(3)
        orbit_choice    = x(4)

        tt = transfer_time * 86400.d0

        ! Initialise velocities

        min_vel = 1.d6

        ! Compute the candidate position

        call CANDIDATE_POSITION(transfer_epoch, state_can)

        ! Main iteration loop: go through the data file and compute Lamberts to that state

        ! Get states via interpolation of the dataset

        call BSPLINE_INTERPOLATE(n_mnfd, orbit_choice, state_targ)

        ! print *, state_targ

        ! Integrate this backwards by t_end

        call integrate(state_targ, t_end)

        ! print *, state_targ

        ! Rotate into the global frame

        call GLOBAL_ROTATE(state_targ, transfer_epoch+tt, state_dim) ! 0 => Rotate in at zero ephemeris seconds (tied to zero epoch)

        ! Dimensionalise

        state_dim(1:3) = state_dim(1:3) * position_dimensionalise_quotient
        state_dim(4:6) = state_dim(4:6) * velocity_dimensionalise_quotient

        state_rot = state_dim

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

    RETURN

END

end program main
