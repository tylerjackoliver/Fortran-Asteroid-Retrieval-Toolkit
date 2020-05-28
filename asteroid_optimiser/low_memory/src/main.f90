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
    use local_optimisation

    implicit none

    integer             :: target_count
    integer, parameter  :: number_of_desired_solutions = 5
    integer             :: number_of_solutions = 0

    character(7)        :: targ_can_temp
    
    logical             :: first_load = .true.

    double precision    :: min_vel
    double precision    :: test(6)

    ! Initialise MPI functionality; read-in datasets

    call MPI_VARIABLE_INIT()

    ! Perform the chunk of our work for the dataset

    do target_count = (mpi_id_world+1), size(targ_can_array), mpi_world_size

        write(targ_can_temp, '(I7)'), targ_can_array(target_count)
        targ_can = trim(targ_can_temp)

        if (first_load) then
            
            call variable_init()
        
        else

            call intermediate_variable_init()

        end if

        call run_global_optim()

        !
        ! Change number of desired solutions if we don't have that many pareto points!
        !

        number_of_solutions = min(PF(1), real(number_of_desired_solutions))

        call initialise_local_optim(number_of_solutions)
        call run_local_optim()

        call intermediate_variable_destruct()

    end do
        
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    call MPI_VARIABLE_DESTRUCT()

contains

    subroutine run_local_optim()

        use iso_fortran_env, only: output_unit

        integer :: solution_number

        do solution_number = 1, number_of_solutions

            call initialise_local_variables(solution_number) ! Initialise bounds etc
            
            do while (optim_stop .eq. 0)

                ! Evaluate objective function and constraints (none)

                call problem_function(f, xopt)

                ! Call MIDACO

                call midaco(1, O, N, NI, M, ME, XOPT, F, G, XL, XU, optim_flag, optim_stop, &
                            param, rw, lrw, iw, liw, pf, lpf, save_to_file, max_eval, &
                            local_time, print_eval)

            end do

            call get_pareto_front() ! With new PF etc
            
            call destroy_local_variables()

            write(*, '(A, I2, A, I2)') "Completed locally optimising condition ", solution_number, " of ", &
                                        number_of_solutions

        end do

    end subroutine run_local_optim

    subroutine run_global_optim()

        do while (optim_stop .eq. 0) 

            ! Evaluate objective function and constraints (none)

            call problem_function(f, xopt)

            ! Call MIDACO

            call midaco(1, O, N, NI, M, ME, XOPT, F, G, XL, XU, optim_flag, optim_stop, &
                        param, rw, lrw, iw, liw, pf, lpf, save_to_file, max_eval, &
                        max_time, print_eval)

        end do

    end subroutine run_global_optim

    subroutine problem_function(F, X)

        ! Output: F, given contraints G and X

        implicit none

        double precision, intent(out)   :: F(2)
        double precision, intent(in)    :: X(5)

        double precision                :: min_vel

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

        use constants                                                               ! Standardises constants
        use state_determination                                                     ! Use the state determination routines (self-written)
        use fortran_astrodynamics_toolkit                                           ! Use the FAT, third-party astrodynamics library
        use problem_parameters                                                      ! Problem set-up parameters
        use compute_spline
        use variable_initialisation

        implicit none

        double precision		                                :: x(5)                         ! Input state vector
        
        double precision                                        :: state_can(6)                 ! Asteroid candidate state
        double precision		                                :: state_targ(6)                ! Un-rotated target state
        double precision                                        :: state_dim(6)
        double precision		                                :: state_rot(6)                 ! Rotated target state

        double precision                                        :: transfer_epoch               ! Epoch of transfer
        double precision		                                :: tt 							! Transfer time
        double precision                                        :: t_end                        ! Desired backwards integration time
        double precision                                        :: n_mnfd                       ! Desired point along the orbit
        double precision                                        :: orbit_choice                 ! The orbit we're using

        double precision		                                :: min_vel                      ! Minimum velocity transfer

        ! Initialise variables from input vector

        transfer_epoch  = x(1)
        tt              = x(2)
        t_end           = x(3)
        n_mnfd          = x(4)
        orbit_choice    = x(5)

        ! if (x(1) .lt. 0 .or. x(2) .lt. 0 .or. x(4) .lt. 0) then

        !     min_vel = 1.d6
        !     RETURN

        ! end if

        ! Compute the candidate position

        call CANDIDATE_POSITION(transfer_epoch, state_can)

        ! Get states via interpolation of the dataset

        call BSPLINE_INTERPOLATE(n_mnfd, orbit_choice, state_targ)

        ! Integrate this backwards by t_end

        call integrate(state_targ, t_end)

        ! Rotate into the global frame

        call GLOBAL_ROTATE(state_targ, transfer_epoch+tt, state_rot) ! Rotate in at correct ephemeris seconds (tied to zero epoch)

        ! Dimensionalise

        state_dim(1:3) = state_rot(1:3) * position_dimensionalise_quotient
        state_dim(4:6) = state_rot(4:6) * velocity_dimensionalise_quotient

        !
        ! Now, call the lambert solver on state_rot (the target
        ! state rotated to the correct epoch). First, with long_way
        ! set to false.
        !

        call lambert(state_can, state_dim, tt, min_vel)

        RETURN

    END

    subroutine lambert(r1, r2, tt, mu, min_vel)

        double precision, intent(in)    :: r1(:)                ! Beginning of Lambert arc
        double precision, intent(in)    :: r2(:)                ! End of Lambert arc
        double precision, intent(in)    :: tt                   ! Length of Lambert arc (related time unit)
        double precision, intent(in)    :: mu                   ! Gravitational parameter of host body

        double precision, intent(out)   :: min_vel              ! Minimum velocity transfer from multi_rev, and long_way flags

        double precision                :: transfer_v1          ! Velocity needed to get onto the Lambert arc
        double precision                :: transfer_v2          ! Velocity needed to get off of the Lambert arc
        double precision                :: transfer_vel         ! Total transfer velocity

        double precision                :: vcx, vcy, vcz        ! Candidate velocity
        double precision                :: vtx, vty, vtz        ! Target velocity
        double precision                :: vx1, vy1, vz1        ! Velocity at start of Lambert
        double precision                :: vx2, vy2, vz2        ! Velocity at end of Lambert

        double precision, allocatable   :: v1(:, :), v2(:, :)   ! Velocity at extremum of Lambert

        integer                         :: num_rows             ! Number of solutions
        integer                         :: multi_rev = 4        ! Number of revolutions to consider in the arc (up to, not including)
        integer                         :: i, j                 ! Iteration sentinels

        logical                         :: run_ok               ! Success flag
        logical                         :: long_way             ! Which way round the Lambert?

        ! Initialise minimum velocity; set trivially high to force update

        min_vel = 1.d6

        long_way = .false.

        ! Define initial and final velocities of the target and of
        ! the candidate

        vcx = r1(4)                                                                                                  ! Candidate; from arguments
        vcy = r1(5)                                                                                                  ! Candidate; from arguments
        vcz = r1(6)                                                                                                  ! Candidate; from arguments

        vtx = r1(4)                                                                                                  ! Rotated state: from file
        vty = r2(5)                                                                                                  ! Rotated state: from file
        vtz = r2(6)                                                                                                  ! Rotated state: from file

        call solve_lambert_izzo(r1(1:3), r2(1:3), tt, mu, long_way, multi_rev, v1, v2, run_ok)

        if (run_ok .eqv. .false.) then

            write(*,*) "Warning: Izzo's solver failed."

            ! We will re-run the Lambert problem using Gooding's
            ! algorithm, which is typically more robust (at the
            ! expense of speed.)

            call solve_lambert_gooding(state_can(1:3),state_dim(1:3),tt,mu,long_way,multi_rev,v1,v2,run_ok)

            if (run_ok .eqv. .false.) then

                min_vel = 1d6

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

        long_way = .true.

        call solve_lambert_izzo(state_can(1:3), state_dim(1:3),tt, mu, long_way,multi_rev,v1,v2,run_ok)

        if (run_ok .eqv. .false.) then

            write(*,*) "Warning: Izzo's solver failed."

            ! We will re-run the Lambert problem using Gooding's
            ! algorithm, which is typically more robust (at the
            ! expense of speed.)

            call solve_lambert_gooding(state_can(1:3),state_dim(1:3),tt,mu,long_way,multi_rev,v1,v2,&
                run_ok)
                
            if (run_ok .eqv. .false.) then

                min_vel = 1d6

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

    end subroutine lambert

end program main
