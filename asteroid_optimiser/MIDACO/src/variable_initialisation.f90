module variable_initialisation

    use precision_kinds
    use constants
    use state_determination
    use problem_parameters

    integer, parameter      :: O = 2                                                        ! Number of objectives
    integer, parameter      :: N = 4                                                        ! Number of variables
    integer, parameter      :: NI = 0                                                       ! Number of integer variables
    integer, parameter      :: M = 0                                                        ! Number of contrains
    integer, parameter      :: ME = 0                                                       ! Number of equality constraints
    integer, parameter      :: LIW = 5000                                                   ! Integer workspace
    integer, parameter      :: LRW = 20000                                                  ! Real workspace
    integer, parameter      :: LPF = 2000000                                                  ! Number of pareto front

    double precision        :: XL(4), XU(4)                                                 ! Upper and lower bounds
    double precision        :: XOPT(4)                                                      ! OPtimisation variables
    double precision        :: F(2)                                                         ! Objectives
    double precision        :: G(4)                                                         ! Constraint arrays; initialised but never used
    double precision        :: PARAM(13)                                                    ! MIDACO parameters
    double precision        :: RW(LRW), PF(LPF)                                             ! Workspace and pareto front

    double precision, allocatable :: pareto_front_data(:, :)                                ! Stores the custom pareto_front_data for generation of ancillary data

    integer                 :: optim_flag                                                   ! Optimiser information flag
    integer                 :: optim_stop                                                   ! Optimiser stopping variable
    integer                 :: IW(LIW)                                                      ! Integer workspace
    integer                 :: max_time                                                     ! Maximum walltime for optimisations
    integer                 :: max_eval                                                     ! Maximum function evaluations
    integer                 :: print_eval                                                   ! How often to print
    integer                 :: save_to_file                                                 ! Output verbosity
    integer                 :: number_of_pareto_solutions                                   ! ...number of Pareto solutions

    character*60            :: key                                                          ! License key bit

    double precision        :: time_lower                                                   ! Lower bound for the optimisation time - dummy here
    double precision        :: time_upper                                                   ! Upper bound for the optimisation time - dummy here

    logical                 :: is_loaded = .false.

    integer                 :: num_targets                                                  ! Number of items in the input file
    integer                 :: num_orbits                                                   ! Number of orbits in the file

    integer, parameter      :: iunit = 70                                                       ! Unit for the Pareto file
    
    integer, parameter      :: t_end_disc = 100                                             ! Parameterisation in backwards time
    integer, parameter      :: n_mnfd_disc = 360                                            ! Parameterisation around the orbit

    contains

        subroutine variable_init()

            call nice_printing()

            ! License key

            key = 'MIDACO_LIMITED_VERSION___[CREATIVE_COMMONS_BY-NC-ND_LICENSE]'

            ! Load SPICE kernels

            call FURNSH('../data/de414.bsp')                    
            call FURNSH('../data/naif0008.tls')
            call FURNSH('../data/'//targ_can//'.bsp')

            ! Initialise transfer epoch bounds - get the state of the candidate at correct epoch

            call GET_STATE(targ_can, time_lower, time_upper)

            ! Optimiser bounds

            XL(1) = time_lower                                                             ! Minimum transfer epoch (ephemeris seconds)
            XL(2) = minimum_transfer_time                                                  ! Minimum transfer duration (seconds)
            XL(3) = 1                                                                      ! t_end
            XL(4) = 1                                                                      ! n_mnfd
            XU(1) = time_upper                                                             ! Maximum transfer epoch (ephemeris seconds)
            XU(2) = maximum_transfer_time                                                  ! Maximum transfer epoch (seconds)
            XU(3) = 100                                                                    ! t_end
            XU(4) = 360                                                                    ! n_mnfd
        
            ! Starting point, XOPT

            XOPT = (/(time_lower+time_upper)*.5d0, 750.d0 * 86400d0, 2.d0, 90.d0/) ! EXACT MIDDLE OF THE SET

            ! Maximum function evaluations

            max_eval = 99999999
            max_time = 60 * 60 * 24 * 1 ! 2 days

            ! Printing options

            print_eval = 100
            save_to_file = 0 ! Save solution to text files

            ! Choose MIDACO parameters (FOR ADVANCED USERS)

            PARAM( 1) = 0.0D0       ! ACCURACY
            PARAM( 2) = 0.0D0       ! SEED
            PARAM( 3) = 0.0D0       ! FSTOP
            PARAM( 4) = 0.0D0       ! ALGOSTOP
            PARAM( 5) = 0.0D0       ! EVALSTOP
            PARAM( 6) = 0.D0       ! FOCUS
            PARAM( 7) = 0.0D0       ! ANTS
            PARAM( 8) = 0.0D0       ! KERNEL
            PARAM( 9) = 0.0D0       ! ORACLE
            PARAM(10) = 100000D0      ! PARETOMAX
            PARAM(11) = 0.000001D0   ! EPSILON  
            PARAM(12) = -1.000D0     ! BALANCE => focus only on first objective function (DeltaV, not tt)
            PARAM(13) = 0.0D0       ! CHARACTER 

            if (.not. is_loaded) call load_data()

        end subroutine variable_init

        subroutine variable_destruct()

            ! ////////////////////////////////////////////////////////////
            !
            ! Ensures graceful program closure
            !
            ! ////////////////////////////////////////////////////////////

            deallocate(dataset)
            deallocate(perturbed_conds_dataset)

            ! Unload SPICE kernels
                    
            call UNLOAD('../data/de414.bsp')                    
            call UNLOAD('../data/naif0008.tls')
            call UNLOAD('../data/'//targ_can//'.bsp')

            close(iunit)

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
            integer     :: i, j, k
            real        :: dum

            if (is_loaded) error STOP "Database has already been loaded. Likely program logic error. Aborting..."

            ! Open file; set number of targets (== number of lines) to zero

            write(*, '(A)', advance="no") " Determining number of orbits in the data-file..."
            num_targets = 0
            open(97, file=datafile)
            open(98, file=original_orbit_data)
            open(iunit, file="../data/paretoFront1_"//targ_can)

            ! Read through in do-loop until EOF to determine number of lines

            do

                read(97,*,iostat=io_state)                                                  ! Not interested in contents yet
                
                if (io_state .lt. 0) then                                                   ! If an error occured (likely EOF)

                    exit                                                                    ! Do loop

                end if

                num_targets = num_targets + 1

            end do

            ! Move from number of targets -> number of orbits

            num_orbits = num_targets / (n_mnfd_disc * t_end_disc)

            write(*, *) "done."
            print *, "Number of target points is ", num_targets
            print *, "Number of orbits is ", num_orbits
            print *, "Number of n_mnfd is ", n_mnfd_disc
            print *, "Number of t_end is ",  t_end_disc

            ! Now we know the number of lines, rewind the file pointer and allocate the database
            ! file into memory (likely very large, so approach with caution on low-memory machines)

            rewind(97)

            allocate(dataset(t_end_disc, n_mnfd_disc, num_orbits, 7))   ! (t_end x n_mnfd x J x 6)
            allocate(perturbed_conds_dataset(n_mnfd_disc, num_orbits, 6))                ! (n_mnfd x J x 7)

            ! Now read the data into memory

            write(*, '(A)', advance="no") " Reading dataset into memory..."

            do k = 1, num_orbits

                do j = 1, n_mnfd_disc

                    do i = 1, t_end_disc

                        read(97, *) dataset(i, j, k, :) ! Dum == time

                    end do ! i

                    read(98, *) perturbed_conds_dataset(j, k, :) ! Dum == time

                end do ! j

            end do ! k

            print *, "done."

            ! Close file; set is_loaded to true.

            close(97)
            close(98)

            is_loaded = .true.

        end subroutine load_data

        subroutine nice_printing()

            print *, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            print *, " "
            print *, "  WELCOME TO THE FORTRAN ASTEROID RETRIEVAL TOOL (FART)  "
            print *, " "
            print *, "                  Author: Jack Tyler                     "
            print *, "              Email: jack.tyler@soton.ac.uk              "
            print *, " "
            print *, " This version last updated 2020-02-03                    "
            print *, " "
            print *, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            print *,
            print *, "Parameters specified"
            print *, "~~~~~~~~~~~~~~~~~~~~"
            print *,
            print *, "Input data-file: ", datafile
            print *, "Original orbit data:", original_orbit_data
            print *, "Target candidate: ", targ_can
            print *,
            print *, "Beginning run "
            print *, "~~~~~~~~~~~~~ "

        end subroutine nice_printing


end module variable_initialisation
        
