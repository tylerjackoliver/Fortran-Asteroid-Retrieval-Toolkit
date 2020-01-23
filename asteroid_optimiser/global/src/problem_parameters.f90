module problem_parameters
    
    use precision_kinds
    use constants
    use state_determination

    implicit none

    ! ////////////////////////////////////////////////////////////
    !
    ! Contains variables that control the program runtime
    !
    ! Inputs
    ! ~~~~~~
    ! targ_can: SPK ID of the target candidate
    ! datafile: Location of target dataset file
    ! dataset: Allocatable double-precision array that contains the target section
    ! is_loaded: Whether the dataset has been loaded into memory
    ! minimum_transfer_time: Minimum duration of Lambert arc (optimiser constraint)
    ! maximum_transfer_time: Maximum duration of Lambert arc (optimiser constraint)
    !
    ! ///////////////////////////////////////////////////////////

    ! //////////////////////////////
    ! 
    ! PROBLEM PARAMETERS - CHANGE ME
    !
    character(*), parameter :: targ_can = '3390109'
    character(*), parameter :: datafile = '../data/2020-01-12_L2PlanarBackCondsGlobal.csv'
    real*8,       parameter :: minimum_transfer_time = 365.d0* 86400                        ! Seconds
    real*8,       parameter :: maximum_transfer_time = 1600.d0 * 86400                      ! Seconds
    !
    ! //////////////////////////////


    real*8,     allocatable :: dataset(:,:) 
    
    real*8                  :: X0(30,20), F0(20), MIN(2), MAX(2), VALUE                     ! Optimiser variables
    real*8                  :: state_can_orig(6)                                            ! Un-rotated asteroid candidate state
    real*8                  :: state_epoch                                                  ! Epoch of the candidate state
    real*8                  :: time_lower                                                   ! Lower bound for the optimisation time - dummy here
    real*8                  :: time_upper                                                   ! Upper bound for the optimisation time - dummy here

    logical                 :: is_loaded = .false.

    integer                 :: M, NPARM, NSIG, NSAMPL, NSEL, IPR                            ! Optimiser variables
    integer                 :: num_targets                                                  ! Number of items in the input file

    contains

        subroutine variable_init()

            ! Number of significant figures in optimiser print-outs
        
            NSIG   = 9
        
            ! Various optimiser variables
        
            M      = 1                                                                      ! No longer used
            NPARM  = 2                                                                      ! Problem dimensions
            NSAMPL = 15                                                                     ! Number of points for sampling. 50 recommended in documentation
            NSEL   = 2                                                                      ! Number of points selected for starting points in sampling
            IPR    = 77                                                                     ! File unit used for writing output

            ! Initialise transfer epoch bounds - get the state of the candidate at correct epoch

            call STATE_FINDER(targ_can, state_can_orig, state_epoch, time_lower, time_upper)

            ! Optimiser bounds

            call FURNSH('../data/naif0008.tls')
            call STR2ET('2032 Feb 15 00:00', time_lower)
            call STR2ET('2032 Mar 15 00:00', time_upper)
            call UNLOAD('../data/naif0008.tls')

            MIN(1) = time_lower                                                             ! Minimum transfer epoch (ephemeris seconds)
            MIN(2) = 1200.D0 * 86400.D0                                                       ! Minimum transfer duration (seconds)
            MAX(1) = time_upper                                                             ! Maximum transfer epoch (ephemeris seconds)
            MAX(2) = 1300.D0 * 86400.D0                                                     ! Maximum transfer epoch (seconds)
        
            if (.not. is_loaded) call load_data()

            ! Open the output file for the optimiser

            OPEN(IPR, FILE='OUTPUT')

        end subroutine variable_init

        subroutine variable_destruct()

            ! ////////////////////////////////////////////////////////////
            !
            ! Ensures graceful program closure
            !
            ! ////////////////////////////////////////////////////////////

            close(IPR)
            deallocate(dataset)

        end subroutine variable_destruct

        subroutine load_data()

            ! ////////////////////////////////////////////////////////////
            ! 
            ! Called if and only if is_loaded is false: loads the datafile
            ! into memory to be accessed by all routines, and sets num_targets
            ! accordingly
            !
            ! ////////////////////////////////////////////////////////////

            integer     :: io_state                                                         ! Tracks EOF
            integer     :: i

            if (is_loaded) then

                error STOP "Database has already been loaded. Likely program logic error. Aborting..."

            end if

            ! Open file; set number of targets (== number of lines) to zero

            write(*, '(A)', advance="no") "Determining number of lines..."

            num_targets = 0
            open(37, file=datafile)

            ! Read through in do-loop until EOF to determine number of lines

            do

                read(37,*,iostat=io_state)                                                  ! Not interested in contents yet
                
                if (io_state .lt. 0) then                                                   ! If an error occured (likely EOF)

                    exit                                                                    ! Do loop

                end if

                num_targets = num_targets + 1

            end do

            write(*, '(A)') "done."

            ! Now we know the number of lines, rewind the file pointer and allocate the database
            ! file into memory (likely very large, so approach with caution on low-memory machines)

            rewind(37)
            allocate(dataset(num_targets, 6))

            ! Now read the data into memory

            write(*, '(A)', advance="no") "Reading dataset into memory..."

            do i = 1, num_targets

                read(37, *) dataset(i, :)

            end do

            write(*, '(A, I8)') "done. Number of targets is", num_targets

            ! Close file; set is_loaded to true.

            close(37)
            is_loaded = .true.

        end subroutine load_data

end module problem_parameters
