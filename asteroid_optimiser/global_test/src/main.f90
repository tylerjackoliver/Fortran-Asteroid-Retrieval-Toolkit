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

    ! Load variables into memory

    call VARIABLE_INIT()

    ! Call the optimiser

    call GLOBAL(MIN, MAX, NPARM, M, NSAMPL, NSEL, IPR, NSIG, X0, NC, F0)

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
        use lambert_transfer_module

        implicit none

        integer, intent(in)                         :: NPARMTR, MM
            
        real*8		                                :: x(:)                         ! Input state vector
        real*8                                      :: state_can(6)                 ! Asteroid candidate state
        real*8		                                :: state_targ(6)                ! Un-rotated target state
        real*8		                                :: state_rot(6)                 ! Rotated target state
        real*8                                      :: transfer_epoch               ! Epoch of transfer
        real*8		                                :: tt 							! Transfer time
        real*8		                                :: transfer_vel                 ! Transfer velocities
        real*8		                                :: min_vel                      ! Minimum velocity transfer

        logical                                     :: long_way						! Which direction for the Lambert transfer
        logical                                     :: run_ok                       ! Boolean success variable for the Lambert

        integer		                                :: multi_rev = 4                ! Number of Lambert arc revolutions (up to)
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

            call ROTATOR(state_targ, transfer_epoch, tt, state_rot)

            call LAMBERT_TRANSFER(state_can, state_rot, tt, multi_rev, transfer_vel)

            if (transfer_vel < min_vel) then

                min_vel = transfer_vel
                best_index = itercount

            end if

        end do main_loop

    rewind(69)
    close(69)

    print *, "Function evaluation completed."

    RETURN

END
