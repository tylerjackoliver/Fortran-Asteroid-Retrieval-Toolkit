! This program implements the Lambert transfer methodology, originally
! written in MATLAB, for Fortran routines. 
!
! This progam is to form one of the main routines for the topimisation
! process in Fortran, and as such a quick performance is paramount.
! Thus, this code will be translated into a parallel implementation at
! some point in the near future.
!
! This program has the following dependencies:
!
!   - Fortran Astrodynamics Toolkit: lambert_solver.f90 module, which
!   provides the Izzo and Gooding solvers. This module has NOT yet been
!   stress-tested against its equivalent PyKep implementation, although
!   this is planned in the very near future.
!
!   - Precision_kinds module: defines the different types of precision
!   for numbers of type REAL and DOUBLE PRECISION, and enforces that
!   double and quadruple precision are indeed twice and quadruple the
!   precision of a single-precision integer, respectively.
!
!   - Constants module: provides useful constants, declared in double
!   precision, for use during the solving routines.
!   
!   - Input data files: This/these files provide the target states for
!   the computation of the Lambert transfers. For the time being, this
!   is being kept as a single large file, but when parallel processing
!   is implemented, this may be subject to change. The location of this
!   file is specified by the CHAR variable 'input_file'.
!
!   - The compiler being used MUST supoprt OpenMP for any of the
!   multi-threading improvements to take place. While the only cluster
!   available is Lyceum2, this program will not use MPI as no further
!   performance increases may be made. This will likely change in the
!   future.
!
!   Brief Variable descriptions
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!   - Transfer_epoch: Time (epehemeris seconds past J2000) at which the
!   transfer occurs
!   - Transfer_time: Length of transfer (seconds)
!

    program main

        use precision_kinds                                       ! Defines appropriate variable kinds for double and quadruple precision
        use constants                                             ! Standardises constants
        use state_determination                                   ! Use the state determination routines (self-written)
        use fortran_astrodynamics_toolkit                         ! Use the FAT, third-party astrodynamics library
    !   use omp_lib

        implicit none

        integer(kind=sp)  :: index                                ! Target state in input file
        real(kind=dp)  :: transfer_time_looper, transfer_epoch_looper

        real(kind=dp)     :: state_can(6)                         ! Asteroid candidate state
        real(kind=dp)     :: state_can_orig(6)                    ! Un-rotated asteroid candidate state
        real(kind=dp)     :: min_vel                              ! Minimum velocity transfer found
        real(kind=dp)     :: transfer_epoch                       ! Epoch of transfer
        real(kind=dp)     :: transfer_time                        ! Time-of-flight
        real(kind=dp)     :: state_epoch                          ! Epoch of the candidate state
        real(kind=dp)     :: time_lower                           ! Lower bound for the time search time
        real(kind=dp)     :: time_upper                           ! Upper bound " "

        real(kind=dp)     :: min_transfer_epoch
        real(kind=dp)     :: min_transfer_time
        real(kind=dp)     :: min_loop_vel

        character(*), parameter  :: targ_can='3435539'                   ! Target candidate string
        character(*), parameter :: input_file='../data/planarL2_newRange_globalBackwards.dat'
        ! character(*), parameter  :: input_file = '../data/states_globalL2Planar.csv'

        ! //////////////////////////////////////////////////////////
        ! Now, use the transfer time and epoch to determine the state
        ! of the candidate at this position.
        !
        ! For now, this is implemented using the SPICE library in the
        ! get_state() subroutine, but this may be subject to change to
        ! a universal variables approach (using PROP2B) if a more
        ! direct replication of the original opimisation process is
        ! required.
        !
        ! EDIT: This now implements PROP2B, a universal variables
        ! approach
        ! ///////////////////////////////////////////////////////////

        call TRANSFER_CALC

    end program main


    subroutine CANDIDATE_POSITION(t_epoch,state_epoch,orig_state,state_can)

        use precision_kinds
        use constants

        REAL(kind=dp), intent(in)  :: t_epoch                               ! Initial time
        REAL(kind=dp), intent(in)  :: state_epoch                           ! Epoch for the desired state
        REAL(kind=dp), intent(in)  :: orig_state(6)                         ! Original state
        
        REAL(kind=dp), intent(out) :: state_can(6)                          ! Final state

        ! Call PROP2B with the propagation time being the difference
        ! between t_epoch and state_epoch

        call PROP2B(mu, orig_state, t_epoch - state_epoch, state_can)

    end subroutine CANDIDATE_POSITION


    subroutine TRANSFER_CALC

        use precision_kinds
        use constants
        use state_determination
        use omp_lib
        use cudafor
        
        REAL(kind=dp)                              :: t                   ! Transfer epoch
        REAL(kind=dp)                              :: tt                  ! Transfer time
        REAL(kind=dp)                              :: min_vel             ! Minimum-velocity transfer found
        REAL(kind=dp)                              :: transfer_time
        REAL(kind=dp)                              :: transfer_epoch 
        REAL(kind=dp)                              :: min_loop_vel

        INTEGER(kind=sp)                           :: f_index             ! State it corresponds to
        INTEGER(kind=sp)                           :: num_iters = 0       ! State it corresponds to
    
        REAL(kind=dp)                              :: opt_state(6)        ! Optimum target state found; debugging
        REAL(kind=dp)                              :: opt_state_unrot(6) 
        REAL(kind=dp)                              :: time_start, time_finish

        real(kind=dp)                              :: state_can(6)        ! Asteroid candidate state
        real(kind=dp)                              :: state_can_orig(6)   ! Un-rotated asteroid candidate state
        real(kind=dp)                              :: state_epoch         ! Epoch of the candidate state
        real(kind=dp)                              :: time_lower          ! Lower bound for the time search time
        real(kind=dp)                              :: time_upper          ! Upper bound " "

        REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: state_targ_arr      ! State of the target (global 3BP)
        REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE :: state_rot_arr       ! Rotated states of the target

        REAL(kind=dp), DIMENSION(:),   ALLOCATABLE :: min_vel_arr         ! Minimum velocity of target

        INTEGER(kind=sp)                           :: io_state            ! IO error variable
        INTEGER(kind=sp)                           :: itercount           ! Iteration counter

        INTEGER(kind=sp)                           :: num_lines           ! Number of lines in the input file
        INTEGER(kind=dp)                           :: i                   ! Loop counter

        REAL(kind=dp)                           :: transfer_epoch_loopercd  
        REAL(kind=dp)                           :: transfer_time_looper

        character(*), parameter  :: targ_can='3435539'                   ! Target candidate string
        character(*), parameter :: input_file='../data/planarL2_newRange_globalBackwards.dat'
        ! character(*), parameter  :: input_file = '../data/states_globalL2Planar.csv'

        ! CUDA variables

        integer                                     :: deviceID
        integer                                     :: numberOfSMs
        integer                                     :: threadsPerBlock
        integer                                     :: cudaError

        ! //////////////////////////////////////////////////////////
        ! Now, use the transfer time and epoch to determine the state
        ! of the candidate at this position.
        !
        ! For now, this is implemented using the SPICE library in the
        ! get_state() subroutine, but this may be subject to change to
        ! a universal variables approach (using PROP2B) if a more
        ! direct replication of the original opimisation process is
        ! required.
        !
        ! EDIT: This now implements PROP2B, a universal variables
        ! approach
        ! ///////////////////////////////////////////////////////////

        call OMP_SET_NUM_THREADS(4)

        transfer_time = 1507.d0*86400.d0

        min_loop_vel = 10.d6

        ! Find the state of the target and the original

        call STATE_FINDER(targ_can,state_can_orig,state_epoch,time_lower,time_upper)

        ! Get the transfer epoch

        call STR2ET('23 Sep 2036 00:00',transfer_epoch)

        transfer_epoch_looper = transfer_epoch 
        transfer_time_looper = transfer_time

        ! Get CUDA information

        cudaError = cudaGetDevice(deviceID)
        cudaError = cudaDeviceGetAttribute(numberOfSMs, 16, deviceID)
        cudaError = cudaDeviceGetAttribute(threadsPerBlock, 1, deviceID)

        ! Open target data file

        open(unit=99,file=input_file,iostat=io_state)

        open(unit=66, file="results.dat")

        write(66, *) "Epoch (s past J2000), Transfer time (s), Velocity (km/s)"

        ! Check for error in opening the file

        if (io_state /= 0) then

            error stop "Error when reading file. Aborting..."

        end if

        ! Print progress message

        write(*,'(A)', advance='no') 'File opened successfully. Determining the number of lines...'

        ! Get the number of lines in the file

        do

            read(99, *, iostat=io_state)

            if (io_state .ne. 0) then

                exit

            end if

            num_lines = num_lines + 1

        end do

        ! Allocate state_targ array

        write(*,*) "Allocating and reading..."

        ! Read the file line-by-line and allocate this into the array
        ! state_targ.
        ! Done in a do-loop until io_state /= 0: at this point, the
        ! loop is terminated.
        !
        ! /// May be quicker to do contiguous read-in

        ! // Contiguous read-in

        allocate(state_targ_arr(num_lines, 6))
        allocate(state_rot_arr(num_lines, 6))
        allocate(min_vel_arr(num_lines))

        rewind(99)

        do i = 1, num_lines

            read(99,*) state_targ_arr(i, 1), state_targ_arr(i, 2), state_targ_arr(i, 3), &
            state_targ_arr(i, 4), state_targ_arr(i, 5), state_targ_arr(i, 6)

        end do

        write(*,*) "Done. Iterating..."

        do t = 0.98 * transfer_epoch, 1.02 * transfer_epoch, 0.01 * transfer_epoch

            do tt = 0.98 * transfer_time, 1.02 * transfer_time, 0.01 * transfer_time

                time_start = omp_get_wtime()

                write(*,*) "Rotating..."

                do i = 1, num_lines

                    call ROTATOR(state_targ_arr(i, :), t, tt, state_rot_arr(i, :))

                end do

                write(*,*) "done."

                call CANDIDATE_POSITION(t,state_epoch,state_can_orig,state_can)

                call velocity_loop<<<numberOfSMs, threadsPerBlock>>>(transfer_time, state_targ_arr, state_rot_arr, min_vel_arr, num_lines)

                cudaError = cudaDeviceSynchronize()

                time_finish = omp_get_wtime()

                write(*,*) "timings: ", time_finish-time_start, "s"

            end do

        end do

        close(99)
        close(66)

    end subroutine TRANSFER_CALC


    subroutine ROTATOR(state_in, epoch, tt, state_out)

        ! ////////////////////////////////////////////////////////////////////////////
        !
        ! Applies relations from Sanchez et. al. to "rotate" a target state by means
        ! of altering its longitude
        !
        ! ///////////////////////////////////////////////////////////////////////////

        use precision_kinds
        use constants

        REAL(kind=dp), intent(in)   :: state_in(6)                                                                              ! Input state; pre-rotation
        REAL(kind=dp), intent(in)   :: epoch                                                                                    ! Epoch of the state
        REAL(kind=dp), intent(in)   :: tt                                                                                       ! Transfer time
        
        REAL(kind=dp), intent(out)  :: state_out(6)                                                                             ! Output state; post-rotation

        REAL(kind=dp)               :: elts(8)                                                                                  ! Elements array

        ! Convert from state -> orbital elements
        ! Epoch does nothing here! So still need to rotate
        ! by (epoch + tt) later on

        call OSCELT(state_in, 0., mu, elts)

        ! Apply relations from Sanchez et. al. The following code
        ! block assumes that the time t is inputted as ephemeris
        ! seconds past J2000 (and thus that it is linked innately to
        ! the epoch of the states used in this work.

        if (abs(elts(4)) < 1.0d-6) then

            elts(4) = 0
            elts(5) = elts(5) + (epoch+tt)*(2.d0*pi)/(365.25d0*86400.d0)                                                                   ! Add on phasing equal to the transfer time                                           

        else

            elts(4) = elts(4) + (epoch+tt)*(2.d0 * pi)/(365.25d0*86400.d0)                                                                   ! Add on phasing equal to the transfer time
            elts(5) = 0

        end if

        ! Convert from orbital elements -> state at the later epoch

        call CONICS(elts, epoch + tt, state_out)

    end subroutine ROTATOR

attributes(global) subroutine velocity_loop(tt, state_targ_arr, state_rot_arr, min_vel_arr, len)

    use precision_kinds
    use constants
    use lambert_module

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: state_targ_arr      ! State of the target (global 3BP)
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN)      :: state_rot_arr       ! Rotated states of the target

    REAL(kind=dp), DIMENSION(:),   ALLOCATABLE, INTENT(INOUT)   :: min_vel_arr         ! Minimum velocity of target

    REAL(kind=dp),                              INTENT(IN)      :: tt

    INTEGER(kind=sp),                           INTENT(IN)      :: len

    REAL(kind=dp)                                               :: state_targ(6)       ! State of the target (3BP)
    REAL(kind=dp)                                               :: state_rot(6)        ! State of the target (inertial frame)
    REAL(kind=dp)                                               :: vx1                 ! x-component of velocity at beginning of Lambert arc
    REAL(kind=dp)                                               :: vy1                 ! y-component of velocity at beginning of Lambert arc
    REAL(kind=dp)                                               :: vz1                 ! z-component of velocity at beginning of Lambert arc
    REAL(kind=dp)                                               :: vx2                 ! x-component of velocity at end of Lambert arc
    REAL(kind=dp)                                               :: vy2                 ! y-component of velocity at end of Lambert arc
    REAL(kind=dp)                                               :: vz2                 ! z-component of velocity at end of Lambert arc
    REAL(kind=dp)                                               :: transfer_v1         ! Transfer velocity required at start of arc
    REAL(kind=dp)                                               :: transfer_v2         ! Transfer velocity required at end of arc
    REAL(kind=dp)                                               :: transfer_vel        ! Total transfer velocity
    REAL(kind=dp)                                               :: vcx                 ! x-velocity of the candidate
    REAL(kind=dp)                                               :: vcy                 ! y-velocity of the candidate
    REAL(kind=dp)                                               :: vcz                 ! z-velocity of the candidate
    REAL(kind=dp)                                               :: vtx                 ! x-velocity of the target
    REAL(kind=dp)                                               :: vty                 ! y-velocity of the target
    REAL(kind=dp)                                               :: vtz                 ! z-velocity of the target

    INTEGER(kind=sp)                                            :: inner_count         ! GPU iteration counter
    INTEGER(kind=sp)                                            :: stride              ! GPU stride count
    INTEGER(kind=sp)                                            :: i
    INTEGER(kind=sp)                                            :: num_rows
    INTEGER(kind=sp)                                            :: multi_rev           ! Number of revolutions to consider in the Lambert arc
    INTEGER(kind=sp)                                            :: j

    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE                  :: v1                  ! ALLOCATABLE attribute required by the Fortran Astrodynamics Toolkit
    REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE                  :: v2                  ! ALLOCATABLE attribute required by the Fortran Astrodynamics Toolkit
    
    LOGICAL                                                     :: long_way            ! Whether to compute the 'long solution' in the Lambert arc
    LOGICAL                                                     :: run_ok              ! Boolean status flag for the Lambert subroutine

    inner_count = threadIdx%x + blockDim%x * blockIdx%x 
    stride = gridDim%x * blockDim%x

    min_vel   = 1.d6
    multi_rev = 4

    do i = inner_count, len, stride

        state_can = state_rot_arr(i, 1:6)
        state_rot = state_targ_arr(i, 1:6)

        min_vel_arr(i) = 1.d6

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

            transfer_v1  = ((vx1-vcx)**2.d0+(vy1-vcy)**2.d0+(vz1-vcz)**2.d0)**.5d0                                                     ! deltaV from the candidate and the beginning of the lambert arc
            transfer_v2  = ((vx2-vtx)**2.d0+(vy2-vty)**2.d0+(vz2-vtz)**2.d0)**.5d0                                                     ! deltaV from the target and the end of the lambert arc
            transfer_vel = transfer_v1 + transfer_v2                                                                        ! Total delta v is the sum of both; both > 0

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
            
            if (transfer_vel < min_vel_arr(i)) then

                min_vel_arr(i) = transfer_vel

            end if


        end do

    end do

end subroutine velocity_loop
