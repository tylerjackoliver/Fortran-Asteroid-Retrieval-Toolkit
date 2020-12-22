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
    use midaco_interface

    implicit none

    integer             :: tt_counter
    integer             :: initial_tt

    character(10)       :: tt_str

    logical             :: first_load = .true.

    double precision    :: initial_xopt(4)

    ! Initialise MPI functionality; read-in datasets

    call MPI_VARIABLE_INIT()
    initial_tt = 1

    ! Open the file with the target candidate solutions in it
    do tt_counter = initial_tt, int( maxval(transfer_time_array) )

        transfer_time = transfer_time_array(tt_counter)

        write(tt_str, '(I4)') tt_counter

        ! If this is the first loop iteration, load the datasets used in the optimisation
        if (first_load) then

                call variable_init()
                first_load = .false.

        ! Otherwise, wipe the current MIDACO state, reload all the required data
        else

                call intermediate_variable_init()

        end if

        ! Run the optimisation
        call run_global_optim()


        ! If we're rank zero, re-call the objective function with the best
        ! solution found in the optimisation to get the minimum cost.
        !
        ! Write this result to a file opened in variable_init()
        if (mpi_id_world .eq. 0) then

            call funct(xopt, pareto_front_minimum)
            write(iunit, *) pareto_front_minimum, tt_counter, xopt

        end if

        ! Set one of the initial guesses on the next iteration to be the optimal
        ! solution here
        initial_xopt = xopt

        ! Destroy variables not needed any more (MIDACO variables etc.)
        call intermediate_variable_destruct()

    end do
    
    ! Wait for all cores to finish
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    ! Finish the MPI runtime
    call MPI_FINALIZE(mpi_err)

contains

    !
    ! @brief Runs the global optimisation routines
    ! All storage is local, except for the current transfer time which is inherited from `main`
    !
    subroutine run_global_optim()

        ! MIDACO-internal parallelism variables
        integer             :: variable_num             ! variable_num \in [0, O]
        integer             :: thread_num               ! Current thread ID
        integer             :: P                        ! Number of threads
        integer             :: status(MPI_STATUS_SIZE)  ! Stores MPI_Status flags

        ! Main results arrays
        double precision    :: XOPT_PAR(1000)           ! Array of parameter sets for each thread
        double precision    :: F_PAR(100)               ! Array of cost function values
        double precision    :: G_PAR(100)               ! Array of constraint violations

        ! Result/design arrays for slaves
        double precision    :: DX(N), DF(N)             ! XOPT, F
        double precision    :: EX(N), EF(N)             ! XOPT, F

        ! For randomly assigning initial conditions
        double precision    :: DIFFERENCE_VECTOR(4)     ! FOR XU - XL
        double precision    :: eta                      ! How far along (0 <= x < 1)

        DIFFERENCE_VECTOR = XU - XL                     ! For setting initial conditions
        P = mpi_world_size                              ! Initialised in MPI_VARIABLE_INIT()

        ! If master, initialise MPI and run
        if (mpi_id_world .eq. 0) then

            ! Copy the starting point XOPT into parallel array XOPT_PAR

            ! First guess should be the previous point
            XOPT_PAR(1:variable_num) = initial_xopt

            ! The rest of the values should be randomly assigned
            do thread_num = 2, P

                do variable_num = 1, N

                    call RANDOM_NUMBER(eta)
                    XOPT_PAR((thread_num-1)*N+variable_num) = XL(variable_num) + &
                        DIFFERENCE_VECTOR(variable_num) * eta

                end do

            end do

            ! Now initialise the optimisation
            do while (optim_stop .eq. 0)

                ! Get the current design variables for each thread
                do thread_num = 2, P

                    do variable_num = 1, N

                        DX(variable_num) = XOPT_PAR((thread_num-1)*N+variable_num)

                    end do

                    ! Send this design variable to the slave to be evaluated
                    CALL MPI_SEND(DX, N, MPI_DOUBLE_PRECISION, thread_num-1, 1, &
                        MPI_COMM_WORLD, mpi_err)

                end do

                ! Evaluate
                call problem_function(F_PAR, XOPT_PAR)

                ! Collect results from slaves
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

        ! Slaves should evaluate the problem function
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
        initial_xopt = XOPT

    end subroutine run_global_optim

    ! @brief Wrap the problem function to the form required in MIDACO
    ! @param[inout] F The cost function values, array of size equal to the number of objectives
    ! @param[inout] X The design parameter values, one-dimensional array of size equal to the number of design variables  
    subroutine problem_function(F, X)

        implicit none

        double precision, intent(out)   :: F(O)
        double precision, intent(in)    :: X(N)

        double precision                :: min_vel

        call funct(x, min_vel)

        F(1) = min_vel

    end subroutine problem_function

    ! @brief Optimiser cost function
    ! @param[inout] X The design parameter values, one-dimensional array of size equal to the number of design variables
    ! @param[inout] MIN_VEL The transfer cost
    SUBROUTINE FUNCT(X, MIN_VEL)

        use constants                                                                           ! Standardises constants
        use state_determination                                                                 ! Use the state determination routines (self-written)
        use fortran_astrodynamics_toolkit                                                       ! Use the FAT, third-party astrodynamics library
        use problem_parameters                                                                  ! Problem set-up parameters
        use compute_spline
        use variable_initialisation
        use integrator

        implicit none

        double precision		                                :: x(N)                         ! Input state vector
        double precision                                        :: state_can(6)                 ! Asteroid candidate state
        double precision		                                :: state_targ(6)                ! Un-rotated target state
        double precision                                        :: state_dim(6)                 ! Dimensionalised state
        double precision		                                :: state_rot(6)                 ! Rotated target state
        double precision                                        :: transfer_epoch               ! Epoch of transfer
        double precision		                                :: tt 							! Transfer time
        double precision                                        :: t_end                        ! Desired backwards integration time
        double precision                                        :: n_mnfd                       ! Desired point along the orbit
        
        double precision                                        :: orbit_choice                 ! The orbit we're using

        double precision		                                :: vx1, vx2, vy1, vy2, vz1, vz2 ! Velocities of the beginning and end of the Lambert arc
        double precision		                                :: vtx, vty, vtz, vcx, vcy, vcz ! Velocities of the candidate and target at beginning and end of Lambert arc
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

        ! Initialise variables from input vector
        transfer_epoch  = x(1)
        t_end           = x(2)
        n_mnfd          = x(3)
        orbit_choice    = x(4)
        tt = transfer_time * 86400.d0                                                           ! Transfer time days => seconds

        ! Initialise velocities
        min_vel = 1.d6

        ! Compute the candidate position - get from ephemeris
        ! Function defined in state_determination.f90
        call CANDIDATE_POSITION(transfer_epoch, state_can)

        ! Get states on the pi/8 plane via interpolation of the dataset
        ! Function defined in compute_spline.f90
        call BSPLINE_INTERPOLATE(n_mnfd, orbit_choice, state_targ)

        ! Integrate this state backwards from the pi / 8 plane by t_end
        call integrate(state_targ, t_end)

        ! Rotate into the global frame
        call GLOBAL_ROTATE(state_targ, transfer_epoch+tt, state_dim) ! Rotate in at the transfer epoch + transfer time in (ephemeris) seconds

        ! Dimensionalise
        ! Constants defined in constants.f90
        state_rot(1:3) = state_dim(1:3) * position_dimensionalise_quotient
        state_rot(4:6) = state_dim(4:6) * velocity_dimensionalise_quotient

        ! Define initial and final velocities of the target and of
        ! the candidate
        vcx = state_can(4)                                                                                                  ! Candidate; from arguments
        vcy = state_can(5)                                                                                                  ! Candidate; from arguments
        vcz = state_can(6)                                                                                                  ! Candidate; from arguments

        vtx = state_rot(4)
        vty = state_rot(5)
        vtz = state_rot(6)

        ! Now, call the lambert solver on state_rot (the target
        ! state rotated to the correct epoch). First, with long_way
        ! set to false.
        long_way = .false.
        call solve_lambert_izzo(state_can(1:3), state_rot(1:3), tt, mu, long_way, multi_rev, v1, v2, run_ok)

        ! If the Izzo Lambert algorithm didn't work, we'll retry with the slower but more robust Gooding solver. 
        if (run_ok .eqv. .false.) then

            write(*,*) "Warning: Izzo's solver failed."

            ! We will re-run the Lambert problem using Gooding's
            ! algorithm, which is typically more robust (at the
            ! expense of speed.)
            call solve_lambert_gooding(state_can(1:3),state_rot(1:3),tt,mu,long_way,multi_rev,v1,v2,run_ok)

            ! If this one failed as well, we've got a problem. Error stop and print an informative message
            if (run_ok .eqv. .false.) then

                error STOP "Lambert arc computation was not successful."

            end if

        end if

        ! Determine the number of solutions; we know the array is num_revs x 3, so get the size here
        num_rows = size(v1, 2)

        ! Go through each solution and compute the capture velocity; compute the velocity for each one and keep the
        ! lowest
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

        ! Determine the number of solutions; we know the array is num_rows x 3, so get the size here
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

    RETURN

END

end program main
