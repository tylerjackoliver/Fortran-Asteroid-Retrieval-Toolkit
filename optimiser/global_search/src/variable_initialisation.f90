! @brief Contains subroutines for initialising, reading and loading data and parallel environments
module variable_initialisation

    use                 mpi
    use                 constants
    use                 state_determination
    use                 problem_parameters
    use, intrinsic   :: iso_c_binding

    implicit none

    integer, parameter                  :: O = 2                                                        ! Number of objectives
    integer, parameter                  :: N = 5                                                        ! Number of variables
    integer, parameter                  :: NI = 0                                                       ! Number of integer variables
    integer, parameter                  :: M = 0                                                        ! Number of contrains
    integer, parameter                  :: ME = 0                                                       ! Number of equality constraints
    integer, parameter                  :: LIW = 5000                                                   ! Integer workspace
    integer, parameter                  :: LRW = 20000                                                  ! Real workspace
    integer, parameter                  :: LPF = 2000000                                                ! Number of pareto front

    double precision                    :: XL(5), XU(5)                                                 ! Upper and lower bounds
    double precision                    :: XOPT(5)                                                      ! OPtimisation variables
    double precision                    :: F(2)                                                         ! Objectives
    double precision                    :: G(1)                                                         ! Constraint arrays; initialised but never used
    double precision                    :: PARAM(13)                                                    ! MIDACO parameters
    double precision                    :: RW(LRW), PF(LPF)                                             ! Workspace and pareto front

    double precision, allocatable       :: pareto_front_data(:, :)                                      ! Stores the custom pareto_front_data for generation of ancillary data

    integer                             :: optim_flag                                                   ! Optimiser information flag
    integer                             :: optim_stop                                                   ! Optimiser stopping variable
    integer                             :: IW(LIW)                                                      ! Integer workspace
    integer                             :: max_eval                                                     ! Maximum function evaluations
    integer                             :: print_eval                                                   ! How often to print
    integer                             :: save_to_file                                                 ! Output verbosity
    integer                             :: number_of_pareto_solutions                                   ! ...number of Pareto solutions
    integer                             :: iostate                                                      ! I/o status flag
    integer                             :: num_candidates                                               ! Number of candidates to search
    integer                             :: candidates_per_processor                                     ! Chunk of candidates per processor
    integer                             :: host                                                         ! Master processor number

    integer, allocatable                :: candidates_to_consider(:)

    character*60                        :: key                                                          ! License key bit

    double precision                    :: time_lower                                                   ! Lower bound for the optimisation time - dummy here
    double precision                    :: time_upper                                                   ! Upper bound for the optimisation time - dummy here

    logical                             :: is_loaded = .false.                                          ! Have we loaded the dataset yet?

    integer(8)                          :: num_targets                                                  ! Number of items in the input file
    integer(8)                          :: num_orbits                                                   ! Number of orbits in the file
    integer(8)                          :: num_orbits_dble                                              ! Long number of orbits in the file
    Integer                             :: iunit                                                        ! Unit for the Pareto file

    ! MPI variables

    integer                             :: mpi_err                                                      ! MPI error integer
    integer                             :: mpi_id_world                                                 ! Core ID of the MPI process
    integer                             :: mpi_id_shared                                                ! ID of the process in the shared-space
    integer                             :: mpi_world_size                                               ! Number of cores in the pool
    integer                             :: mpi_name_string_size                                         ! Length of the name of this processor
    integer                             :: mpi_color                                                    ! For grouping processors by shared memory space

    ! MPI communicators

    integer                             :: node_communicator                                            ! Intra-node communicator
    integer                             :: node_master_communicator                                     ! Master nodes on each core

    ! MPI node naming

    character(MPI_MAX_PROCESSOR_NAME)   :: proc_name                                                    ! Processor name

    ! MPI window variables

    TYPE(C_PTR)                         :: dataset_pointer                                              ! C-style pointer to the shared memory array
    TYPE(C_PTR)                         :: perturbed_conds_pointer                                      ! C-style pointer to the perturbed conds array
    integer                             :: dataset_window                                               ! Window that MPI opens to the shared memory
    integer                             :: perturbed_conds_window                                       ! Window that MPI opens to the shared memory
    integer                             :: info                                                         ! Information integer; MPI
    integer(kind=MPI_ADDRESS_KIND)      :: dataset_size_bytes                                           ! Dataset size in bytes
    integer(8)                          :: row_length                                                   ! Length of a row of the dataset
    integer(8)                          :: double_size                                                  ! Size of a double
    integer                             :: disp_unit                                                    ! Displacement unit for MPI sends
    integer(8)                          :: n_mnfd_disc_dble = n_mnfd_disc                               ! Double-precision integer representation of n_mnfd
    integer(8)                          :: to_send_size                                                 ! Size of message to send
    real*8                              :: dataset_double                                               ! Dataset size in double-precision

    contains

        ! @brief Initialises all the variables needed to open the MPI window
        subroutine mpi_variable_init()

            integer                     :: j, k

            disp_unit=1

            ! Initialise MPI
            call MPI_INIT(mpi_err)
            call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_id_world, mpi_err)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_world_size, mpi_err)
            call MPI_GET_PROCESSOR_NAME(proc_name, mpi_name_string_size, mpi_err)

            if (mpi_id_world .eq. 0) then
                ! Print status message 
                call nice_printing()
                ! Determine the number of orbits in the data file by iterating through once and dividing by n_mnfd_disc
                call determine_num_orbits()
            end if

            ! First, broadcast the number of orbits
            call MPI_BCAST(num_orbits, STORAGE_SIZE(num_orbits), MPI_BYTE, 0, MPI_COMM_WORLD, mpi_err)
            num_orbits_dble = num_orbits

            ! Split the global communicator into smaller communicators,
            ! each of which defines a shared-memory region (each node)
            call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_communicator, mpi_err)
            call MPI_COMM_RANK(node_communicator, mpi_id_shared, mpi_err)
            call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_id_world, mpi_err)

            ! Now create another communicator for *only* the processors that are masters
            ! on their relative nodes
            if (mpi_id_shared .eq. 0) then

               mpi_color = 1               ! Want this to be in the communicator

            else

               mpi_color = MPI_UNDEFINED   ! Do not want this to be in the communicator

            end if

            ! Split the workers into the new communicator
            call MPI_COMM_SPLIT(MPI_COMM_WORLD, mpi_color, 0, node_master_communicator, mpi_err)

            !
            ! Now we have three communicators: WORLD, intra-node, and node-masters
            ! We need to have every node-master have a copy of the target data, so broadcast
            ! the target array *solely* to the node-masters
            !
            ! Allocate memory only for the node masters but NOT the one that has already allocated data
            row_length = 7
            double_size = sizeof(dataset_double)
            
            if (mpi_id_shared .eq. 0) then

                dataset_size_bytes = int(num_orbits * n_mnfd_disc * row_length * double_size, MPI_ADDRESS_KIND)

            else

                dataset_size_bytes = 0

            end if

            ! Allocate a shared memory region on each of the node masters
            call MPI_WIN_ALLOCATE_SHARED(dataset_size_bytes, disp_unit, MPI_INFO_NULL, node_communicator, dataset_pointer, &
                                         dataset_window, mpi_err)
            !
            ! If we're not the master thread, then we won't know the location of the memory space; thus, query the location of the memory here for use later
            ! This returns a C-style pointer
            if (mpi_id_shared /= 0) then

                call MPI_WIN_SHARED_QUERY(dataset_window, 0, dataset_size_bytes, disp_unit, dataset_pointer, mpi_err)

            end if

            ! Now convert the C pointer to a Fortran pointer (syntax for accessing elements is the same as standard ALLOCATEd arrays)
            call C_F_POINTER(dataset_pointer, dataset, (/num_orbits_dble, &
                                                        n_mnfd_disc_dble, &
                                                        row_length/))
  
            ! Load in the data
            ! To avoid racing on the IO, have the global master _only_ read in the data
            if (mpi_id_world .eq. 0) then

                call load_data()
                call load_targets()

            end if

            ! Wait for this data to finish being loaded
            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

            !
            ! Now the global master should BCAST all of the data it just read-in to other node-masters
            ! However, the data it is sending overflows the internal C-type `int` used in MPI, so we can only
            ! send up to 2**32-1 bits at a time.
            !
            ! Here, we eat the overheads and send 'stripes' of the dataset and perturbed_conds dataset out to
            ! the other node-masters. This will not overflow *our* dataset, and makes reassembling the data on the other end
            ! far easier compared to using e.g. MPI structures
            !
            ! Provided all workers call this subroutine, then it is also blocking/synchronous
            ! 
            if (mpi_id_shared .eq. 0) then

                host = 0
                do j = 1, n_mnfd_disc
                    do k = 1, int(num_orbits)
                        ! Broadcast this stripe of the dataset, sending bytes to support arbitrary Fortran position
                        call MPI_BCAST(dataset(k, j, :), STORAGE_SIZE(dataset(k, j, :)/8), MPI_byte, host, &
                                        node_master_communicator, mpi_err)  
                    
                    end do

                end do

            end if
                
            !
            ! Tell all nodes to wait for the node masters to finish, just in case
            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

            !
            ! Now send the number of candidates that we have to all of the other nodes
            ! Master nodes have already read this data in `load_targets()`
            host = 0
            call MPI_BCAST(num_candidates, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)    
            
            if (mpi_id_world .ne. 0) allocate(targ_can_array(num_candidates))
            
            ! Now broadcast the target arrays out
            to_send_size = size(targ_can_array)
            call MPI_BCAST(targ_can_array, size(targ_can_array), MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

        end subroutine mpi_variable_init

        subroutine determine_num_orbits()

            integer :: io_state

            write(*, '(A)', advance="no") "Determining number of lines..."
            num_targets = 0
            open(97, file=datafile)

            ! Read through in do-loop until EOF to determine number of lines
            do

              read(97,*,iostat=io_state)                                                  ! Not interested in contents yet
               
              if (io_state .lt. 0) then                                                   ! If an error occured (likely EOF)

                  exit                                                                    ! Do loop

              end if

              num_targets = num_targets + 1

	        end do

           write(*, '(A)') "done."

           rewind(97)

            num_orbits = num_targets / (n_mnfd_disc)

        end subroutine determine_num_orbits

        subroutine variable_init()

            ! Load SPICE kernels

            character(3) :: run_counter_str

            write(run_counter_str, '(I2)') run_counter

            iunit = 151+mpi_id_world
            call FURNSH('../data/de414.bsp')                    
            call FURNSH('../data/naif0008.tls')
            open(iunit, file="../"//targ_can)

            call STR2ET('Jan 1 2025 00:00', time_lower)
            call STR2ET('Dec 30 2099 00:00', time_upper)

            ! Optimiser bounds

            XL(1) = time_lower                                                              ! Minimum transfer epoch (ephemeris seconds)
            XL(2) = 0.d0                                                                    ! Transfer time
            XL(3) = -25.d0                                                                  ! t_end - 1500 days
            XL(4) = 1                                                                       ! n_mnfd
            XL(5) = 1                                                                       ! J
            XU(1) = time_upper                                                              ! Maximum transfer epoch (ephemeris seconds)
            XU(2) = 1500 * 86400.d0                                                         ! Transfer time
            XU(3) = 0.d0                                                                    ! t_end - at the plane
            XU(4) = n_mnfd_disc-1                                                           ! n_mnfd
            XU(5) = num_orbits-1                                                            ! J

            ! Starting point, XOPT
            XOPT = (XL + XU)/2.d0

            max_eval = 100000

            ! Printing options
            print_eval = 10000
            save_to_file = 0 ! Do not have MIDACO save solution to text files

            ! Choose MIDACO parameters (FOR ADVANCED USERS)

            PARAM( 1) = 0.0D0           ! ACCURACY
            PARAM( 2) = 34.0D0          ! SEED
            PARAM( 3) = 0.0D0           ! FSTOP
            PARAM( 4) = 0.0D0           ! ALGOSTOP
            PARAM( 5) = 0.0D0           ! EVALSTOP
            PARAM( 6) = 0.D0            ! FOCUS
            PARAM( 7) = 500.0D0         ! ANTS
            PARAM( 8) = 100.0D0         ! KERNEL
            PARAM( 9) = 0.0D0           ! ORACLE
            PARAM(10) = 10000D0         ! PARETOMAX
            PARAM(11) = 0.000001D0      ! EPSILON  
            PARAM(12) = -1.000D0        ! BALANCE => focus only on first objective function (DeltaV, not tt)
            PARAM(13) = 0.0D0           ! CHARACTER 

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

        end subroutine variable_init

        ! @brief Remove variables between MIDACO runs
        subroutine intermediate_variable_destruct()

            call UNLOAD(ephemeris_prefix//targ_can//'.bsp')
            close(iunit)

        end subroutine intermediate_variable_destruct

        ! @brief Initialise variables between MIDACO runs
        subroutine intermediate_variable_init()

            use midaco_interface
            implicit none

            integer         :: temp        

            ! File output unit
            iunit = 151+mpi_id_world

            ! Open a file named as the transfer time in days
            open(iunit, file="../"//targ_can)

            call FURNSH(ephemeris_prefix//targ_can//'.bsp')

            ! Optimiser bounds
            call STR2ET('Jan 1 2025 00:00', time_lower)
            call STR2ET('Dec 30 2099 00:00', time_upper)

            XL(1) = time_lower                                                              ! Minimum transfer epoch (ephemeris seconds)
            XL(2) = 0.d0                                                                    ! Transfer time
            XL(3) = -25.d0                                                                  ! t_end - 1500 days
            XL(4) = 2                                                                       ! n_mnfd
            XL(5) = 2                                                                       ! J
            XU(1) = time_upper                                                              ! Maximum transfer epoch (ephemeris seconds)
            XU(2) = 1500 * 86400.d0                                                         ! Transfer time
            XU(3) = 0.d0                                                                    ! t_end - at the plane
            XU(4) = n_mnfd_disc-1                                                           ! n_mnfd
            XU(5) = num_orbits-1                                                            ! J

            ! Starting point, XOPT
            XOPT = (XL + XU)/2.d0

            ! Maximum function evaluations
            max_eval = 99999999
            ! Reset internal MIDACO variables
            optim_flag = 0
            optim_stop = 0

            ! Printing options
            print_eval = 10000
            save_to_file = 0 ! Do not have MIDACO save solution to text files

            ! Choose MIDACO parameters (FOR ADVANCED USERS)

            PARAM( 1) = 0.0D0                       ! ACCURACY
            PARAM( 2) = 34.d0                       ! SEED
            PARAM( 3) = 0.0D0                       ! FSTOP
            PARAM( 4) = 0.0D0                       ! ALGOSTOP
            PARAM( 5) = 0.0D0                       ! EVALSTOP
            PARAM( 6) = 0.D0                        ! FOCUS
            PARAM( 7) = 500.0D0                     ! ANTS
            PARAM( 8) = 100.0D0                     ! KERNEL
            PARAM( 9) = 0.0D0                       ! ORACLE
            PARAM(10) = 100000D0                    ! PARETOMAX
            PARAM(11) = 0.000001D0                  ! EPSILON  
            PARAM(12) = -1.000D0                    ! BALANCE => focus only on first objective function (DeltaV, not tt)
            PARAM(13) = 0.0D0                       ! CHARACTER 

            ! Reset the work arrays and Pareto front arrays
            iw = 0
            rw = 0.d0
            pf = 0.d0

            temp = reset_midaco_run()

        end subroutine intermediate_variable_init

        ! @brief Destroy all MPI variables (pointers, shared memory windows) at the end of program execution
        subroutine mpi_variable_destruct()

            ! ////////////////////////////////////////////////////////////
            !
            ! Ensures graceful program closure
            !
            ! ////////////////////////////////////////////////////////////

            ! Unload SPICE kernels
                    
            call UNLOAD('../data/de414.bsp')                    
            call UNLOAD('../data/naif0008.tls')
            call UNLOAD(ephemeris_prefix//targ_can//'.bsp')

            ! Free the shared-memory space

            call MPI_WIN_FREE(dataset_window, mpi_err)
            call MPI_WIN_FREE(perturbed_conds_window, mpi_err)

            call MPI_FINALIZE(mpi_err)

            close(iunit)

        end subroutine mpi_variable_destruct

        ! @brief Destroy all variables created during the main initialisation process
        subroutine variable_destruct()

            ! ////////////////////////////////////////////////////////////
            !
            ! Ensures graceful program closure
            !
            ! ////////////////////////////////////////////////////////////

            deallocate(dataset)
            ! Unload SPICE kernels
                    
            call UNLOAD('../data/de414.bsp')                    
            call UNLOAD('../data/naif0008.tls')
            call UNLOAD(ephemeris_prefix//targ_can//'.bsp')

            close(iunit)

        end subroutine variable_destruct

        ! @brief Load the reference datasets
        subroutine load_data()

            ! ////////////////////////////////////////////////////////////
            ! 
            ! Called if and only if is_loaded is false: loads the datafile
            ! into memory to be accessed by all routines, and sets num_targets
            ! accordingly
            !
            ! ////////////////////////////////////////////////////////////

            integer     :: j, k

            if (is_loaded) error STOP "Database has already been loaded. Likely program logic error. Aborting..."

            ! Open file; set number of targets (== number of lines) to zero

            write(*, '(A)', advance="no") " Determining number of orbits in the data-file..."
            num_targets = 0

            ! ! Move from number of targets -> number of orbits

            print *, "Number of target points is ", num_targets
            print *, "Number of orbits is ", num_orbits
            print *, "Number of n_mnfd is ", n_mnfd_disc
            rewind(97)

            ! Now read the data into memory

            write(*, '(A)', advance="no") " Reading dataset into memory..."

            do k = 1, int(num_orbits)

                do j = 1, n_mnfd_disc

                        read(97, *) dataset(k, j, :) ! Dum == time

                end do ! j

            end do ! k

            ! Close file; set is_loaded to true.
            close(97)
            is_loaded = .true.

        end subroutine load_data

        ! @brief Load the SPK IDs of the target asteroids
        subroutine load_targets()

            integer                 :: file_unit = 95
            integer                 :: io_state
            integer                 :: candidate

            open(unit=file_unit, file=target_file)

            num_candidates = 0

            do

                read(file_unit, *, iostat=io_state)

                if (io_state .ne. 0) exit

                num_candidates = num_candidates + 1

            end do

            rewind(file_unit)

            allocate(targ_can_array(num_candidates))

            do candidate = 1, num_candidates

                read(file_unit, *) targ_can_array(candidate)

            end do

        end subroutine load_targets


        subroutine nice_printing()

            print *, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            print *, " "
            print *, "  WELCOME TO THE FORTRAN ASTEROID RETRIEVAL TOOL (FART)  "
            print *, " "
            print *, "                  Author: Jack Tyler                     "
            print *, "              Email: jack.tyler@soton.ac.uk              "
            print *, " "
            print *, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            print *, " "
            print *, "Parameters specified"
            print *, "~~~~~~~~~~~~~~~~~~~~"
            print *, " "
            print *, "Input data-file: ", datafile
            print *, "Target candidate: ", targ_can
            print *, " "
            print *, "Beginning run "
            print *, "~~~~~~~~~~~~~ "

        end subroutine nice_printing


end module variable_initialisation
        
