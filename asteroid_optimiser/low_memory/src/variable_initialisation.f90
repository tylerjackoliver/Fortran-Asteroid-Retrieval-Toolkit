!------------------------------------------------------------------------------
! Fortran Asteroid Retrieval Tool (FART) v1.0: variable_initialisation module
!------------------------------------------------------------------------------
!
! MODULE: Variable initialisation
!
!> @author
!> Jack Tyler, University of Southampton
!
! DESCRIPTION: 
!> Contains routines to create parallel working regions based on the MPI shared
!! memory spaces paradig, and distribute and synchronise data across all workers.
!
! REVISION HISTORY:
! 01 Oct 2019 - Initial Version
! 01 Jul 2020 - Refactoring; add Doxygen support
!------------------------------------------------------------------------------

module variable_initialisation

    use                 mpi
    use                 constants
    use                 state_determination
    use                 problem_parameters
    use, intrinsic   :: iso_c_binding

    implicit none

    ! --------------------------------------------
    !
    ! The following variables are unique to MIDACO.
    ! The user should supply their own variables
    ! for their own optimiser where (MIDACO) is given
    ! in the description of the variable.
    !
    ! ---------------------------------------------

    integer, parameter                  :: O = 2                                                        ! Number of objectives
    integer, parameter                  :: N = 5                                                        ! Number of variables
    integer, parameter                  :: NI = 0                                                       ! Number of integer variables
    integer, parameter                  :: M = 0                                                        ! Number of contraints
    integer, parameter                  :: ME = 0                                                       ! Number of equality constraints
    integer, parameter                  :: LIW = 5000                                                   ! Integer workspace (MIDACO)
    integer, parameter                  :: LRW = 20000                                                  ! Real workspace (MIDACO)
    integer, parameter                  :: LPF = 2000000                                                ! Number of pareto front (MIDACO)
    double precision                    :: PARAM(13)                                                    ! MIDACO parameters (MIDACO)
    double precision                    :: RW(LRW), PF(LPF)                                             ! Workspace and pareto front (MIDACO)

    double precision                    :: XL(5), XU(5)                                                 ! Upper and lower bounds
    double precision                    :: XOPT(5)                                                      ! Optimisation variables
    double precision                    :: F(2)                                                         ! Objectives
    double precision                    :: G(1)                                                         ! Constraint arrays; initialised but never used

    double precision, allocatable       :: pareto_front_data(:, :)                                      ! Stores the custom pareto_front_data for generation of ancillary data

    integer                             :: optim_flag                                                   ! Optimiser information flag (MIDACO)
    integer                             :: optim_stop                                                   ! Optimiser stopping variable (MIDACO)
    integer                             :: IW(LIW)                                                      ! Integer workspace (MIDACO)
    integer                             :: max_eval                                                     ! Maximum function evaluations (MIDACO)
    integer                             :: print_eval                                                   ! How often to print (MIDACO)
    integer                             :: save_to_file                                                 ! Output verbosity (MIDACO)
    integer                             :: number_of_pareto_solutions                                   ! ...number of Pareto solutions (MIDACO)

    integer                             :: iostate                                                      ! File I/O status variable
    integer                             :: num_candidates                                               ! Number of candidates considered in the work
    integer                             :: candidates_per_processor                                     ! Number of candidates per processor
    integer                             :: host                                                         ! Host processor ID

    integer, allocatable                :: candidates_to_consider(:)                                    ! Candidates to optimise

    character*60                        :: key                                                          ! License key (MIDACO)

    double precision                    :: time_lower                                                   ! Lower bound for the optimisation time - dummy here
    double precision                    :: time_upper                                                   ! Upper bound for the optimisation time - dummy here

    logical                             :: is_loaded = .false.                                          ! Has the dataset been loaded?

    !
    ! The odd mix of single- and double-precision integers in the variable declarations below is as a result of some variables (particularly the number of lines)
    ! needing to be doubles to prevent overflows. The way MPI is coded for its Fortran/FORTRAN interface means that if one argument is double precision, they all
    ! have to be. How convenient!
    !

    integer(8)                          :: num_targets                                                  ! Number of items in the input file
    integer(8)                          :: num_orbits                                                   ! Number of orbits in the file
    integer(8)                          :: num_orbits_dble

    integer                             :: iunit                                                        ! Unit for the Pareto file

    ! MPI variables

    integer                             :: mpi_err                                                      ! MPI error integer
    integer                             :: mpi_id_world                                                 ! Core ID of the MPI process
    integer                             :: mpi_id_shared                                                ! ID of the process in the shared-space
    integer                             :: mpi_world_size                                               ! Number of workers
    integer                             :: mpi_name_string_size                                         ! Maximum length of MPI name string
    integer                             :: mpi_color                                                    ! MPI inclusion flag in topology set-up

    ! MPI communicators

    integer                             :: node_communicator                                            ! Intra-node communicator
    integer                             :: node_master_communicator                                     ! Master nodes on each core

    ! MPI node naming

    character(MPI_MAX_PROCESSOR_NAME)   :: proc_name                                                    ! Processor name to MPI

    ! MPI window variables

    TYPE(C_PTR)                         :: dataset_pointer                                              ! C-style pointer to the shared memory array
    TYPE(C_PTR)                         :: perturbed_conds_pointer                                      ! C-style pointer to the perturbed conds array
    integer                             :: dataset_window                                               ! Window that MPI opens to the shared memory
    integer                             :: perturbed_conds_window                                       ! Window that MPI opens to the shared memory

    integer                             :: info                                                         ! Another MPI status flag

    integer(kind=MPI_ADDRESS_KIND)      :: dataset_size_bytes                                           ! Size of the target dataset in bytes
    integer(kind=MPI_ADDRESS_KIND)      :: perturbed_conds_size_bytes                                   ! Size of perturbed initial conditions in bytes
    
    integer(8)                          :: row_length = 7                                               ! Length of one row of the target dataset
    integer(8)                          :: row_length_perturbed = 6                                     ! Length of one row of the perturbed conditions dataset
    integer(8)                          :: double_size                                                  ! sizeof(double) holder
    integer                             :: disp_unit                                                    ! Displacement stride, MPI
    integer(8)                          :: n_mnfd_disc_dble = n_mnfd_disc                               ! Number of in-orbit discretisations
    integer(8)                          :: to_send_size                                                 ! Size of message to send
    
    double precision                    :: dataset_double                                               ! 
    double precision                    :: perturbed_conds_double

    contains
    
        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Initialises the MPI topology for the problem, and reads-in and distributes
        !! all the datasets to be used in the problem.
        !
        ! REVISION HISTORY:
        ! 01 Oct 2019 - Initial version
        ! 4 July 2020 - Refactor; add Doxygen support   
        !
        !--------------------------------------------------------------------------- 

        subroutine mpi_variable_init()

            integer                     :: i, j, k

            disp_unit=1

            !
            ! Initialise MPI
            !

            call MPI_INIT(mpi_err)
            call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_id_world, mpi_err)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_world_size, mpi_err)
            call MPI_GET_PROCESSOR_NAME(proc_name, mpi_name_string_size, mpi_err)

            !
            ! If master core, print status message, and determine the number of orbits in the data-file
            !

            if (mpi_id_world .eq. 0) then

                call nice_printing()
                call determine_num_orbits()

            end if

            !
            ! First, broadcast the number of orbits to all other workers
            !

            call MPI_BCAST(num_orbits, STORAGE_SIZE(num_orbits), MPI_BYTE, 0, MPI_COMM_WORLD, mpi_err)

            num_orbits_dble = num_orbits

            !
            ! Split the global communicator into smaller communicators,
            ! each of which defines a shared-memory region (each node)
            !

            call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_communicator, mpi_err)
            call MPI_COMM_RANK(node_communicator, mpi_id_shared, mpi_err)
            call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_id_world, mpi_err)

            !
            ! Now create another communicator for *only* the workers that are masters
            ! on their relative nodes
            !

            if (mpi_id_shared .eq. 0) then

               mpi_color = 1               ! Want this to be in the communicator

            else

               mpi_color = MPI_UNDEFINED   ! Do not want this to be in the communicator

            end if

            call MPI_COMM_SPLIT(MPI_COMM_WORLD, mpi_color, 0, node_master_communicator, mpi_err)

            !
            ! Now we have three communicators: WORLD, intra-node, and node-masters
            ! We need to have every node-master have a copy of the target data, so broadcast
            ! the target array *solely* to the node-masters
            !

            ! Allocate memory only for the node masters but NOT the one that has already allocated data

            double_size = sizeof(dataset_double)
            
            if (mpi_id_shared .eq. 0) then

                dataset_size_bytes = int(num_orbits * n_mnfd_disc * row_length * double_size, MPI_ADDRESS_KIND)
                perturbed_conds_size_bytes = int(num_orbits * row_length_perturbed * n_mnfd_disc * double_size, MPI_ADDRESS_KIND)
            

            else

                dataset_size_bytes = 0
                perturbed_conds_size_bytes = 0

            end if

            call MPI_WIN_ALLOCATE_SHARED(dataset_size_bytes, disp_unit, MPI_INFO_NULL, node_communicator, dataset_pointer, &
                                         dataset_window, mpi_err)
            call MPI_WIN_ALLOCATE_SHARED(perturbed_conds_size_bytes, disp_unit, MPI_INFO_NULL, node_communicator, &
                                         perturbed_conds_pointer, perturbed_conds_window, mpi_err)

            !
            ! If we're not the master thread, then we won't know the location of the memory space; thus, query the location of the memory here for use later
            !

            if (mpi_id_shared /= 0) then

                call MPI_WIN_SHARED_QUERY(dataset_window, 0, dataset_size_bytes, disp_unit, dataset_pointer, mpi_err)
                call MPI_WIN_SHARED_QUERY(perturbed_conds_window, 0, perturbed_conds_size_bytes, disp_unit, &
                                          perturbed_conds_pointer, mpi_err)

            end if

            !
            ! Now convert the C pointer to a Fortran pointer (syntax for accessing elements is the same as standard ALLOCATEd arrays)
            ! 

            call C_F_POINTER(dataset_pointer, dataset, (/num_orbits_dble, &
                                                        n_mnfd_disc_dble, &
                                                        row_length/))
            call C_F_POINTER(perturbed_conds_pointer, perturbed_conds_dataset, &
                (/n_mnfd_disc_dble, num_orbits_dble, row_length_perturbed/))
  
            !
            ! Read in the dataset
            ! To avoid racing on the IO, have the global master _only_ read in the data
            !

            if (mpi_id_world .eq. 0) then

                call load_data()
                call load_targets()

            end if

            ! Wait until the read-in has finished

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

            !
            ! Now the global master should BCAST all of the data it just read-in to other node-masters
            ! However, the data it is sending overflows the internal C-type `int` used in MPI, so we can only
            ! send up to 2**32-1 bits at a time.
            !
            ! Here, we eat the overheads and send 'stripes' of the dataset and perturbed_conds dataset out to
            ! the other node-masters. This will not overflow *our* dataset, and makes reassembling the data on the other end
            ! far easier compared to using e.g. MPI structures (and to be honest, this works even though it's inelegant)
            !
            ! Provided all workers call this subroutine, then it is also blocking/synchronous
            ! 

            if (mpi_id_shared .eq. 0) then

                host = 0

                do j = 1, n_mnfd_disc

                    do k = 1, num_orbits

                        call MPI_BCAST(dataset(k, j, :), STORAGE_SIZE(dataset(k, j, :)/8), MPI_byte, host, &
                                        node_master_communicator, mpi_err)  
                    
                    end do

                end do

                  
                do i = 1, n_mnfd_disc

                        do j = 1, num_orbits

                                ! Send perturbed_conds

                                to_send_size = size(perturbed_conds_dataset(i, j, :))
                                call MPI_BCAST(perturbed_conds_dataset(i,j,:), STORAGE_SIZE(perturbed_conds_dataset(i,j,:)/8), &
                                                MPI_byte, host, node_master_communicator, mpi_err) 

                        end do

                end do

            end if
                
            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

            !
            ! Now send the number of candidates that we have to all of the other nodes
            !

            call MPI_BCAST(num_candidates, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)    
            
            if (mpi_id_world .ne. 0) allocate(targ_can_array(num_candidates))

            to_send_size = size(targ_can_array)
            call MPI_BCAST(targ_can_array, size(targ_can_array), MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

        end subroutine mpi_variable_init

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Determines the number of lines in the input data file.
        !
        ! REVISION HISTORY:
        ! 20 Mar 2018 - Initial version
        ! 3 Jul 2020 - Refactor; Doxygen Support
        !
        !--------------------------------------------------------------------------- 

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

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Initialises optimisation variables, opens target datasets, and loads
        !! ephemerides
        !
        ! REVISION HISTORY:
        ! 10 Oct 2019 - Initial version
        ! 03 Jul 2020 - Refactor; add Doxygen
        !
        !--------------------------------------------------------------------------- 

        subroutine variable_init()

            ! Load SPICE kernels

            iunit = 151+mpi_id_world
            call FURNSH('../data/de414.bsp')                    
            call FURNSH('../data/naif0008.tls')
            call FURNSH(ephemeris_prefix//targ_can//'.bsp')

            !
            ! Open input file
            !

            open(iunit, file="../data/paretoFront_"//targ_can)

            ! Initialise transfer epoch bounds - get the state of the candidate at correct epoch

            call STR2ET('Jan 1 2025 00:00', time_lower)
            call STR2ET('Dec 30 2099 00:00', time_upper)

            ! Optimiser bounds

            XL(1) = time_lower                                                             ! Minimum transfer epoch (ephemeris seconds)
            XL(2) = minimum_transfer_time                                                  ! Minimum transfer duration (seconds)
            XL(3) = -25.d0                                                                  ! t_end - 1500 days
            XL(4) = 1                                                                      ! n_mnfd
            XL(5) = 2                                                                      ! J
            XU(1) = time_upper                                                             ! Maximum transfer epoch (ephemeris seconds)
            XU(2) = maximum_transfer_time                                                  ! Maximum transfer epoch (seconds)
            XU(3) = 0.d0                                                                   ! t_end - at the plane
            XU(4) = n_mnfd_disc                                                            ! n_mnfd
            XU(5) = num_orbits-1                                                           ! J

            ! Starting point, XOPT

            XOPT = (XL + XU)/2.d0

            ! Maximum function evaluations

            max_eval = 100000

            ! Printing options

            print_eval = 1000       ! How often to print (MIDACO)
            save_to_file = 0        ! Save solution to text files (MIDACO)

            ! Choose MIDACO parameters (FOR ADVANCED USERS)

            PARAM( 1) = 0.0D0           ! ACCURACY (MIDACO)
            PARAM( 2) = 0.0D0           ! SEED (MIDACO)
            PARAM( 3) = 0.0D0           ! FSTOP (MIDACO)
            PARAM( 4) = 0.0D0           ! ALGOSTOP (MIDACO)
            PARAM( 5) = 0.0D0           ! EVALSTOP (MIDACO)
            PARAM( 6) = 100.D0          ! FOCUS (MIDACO)
            PARAM( 7) = 100.0D0         ! ANTS (MIDACO)
            PARAM( 8) = 50.0D0          ! KERNEL (MIDACO)
            PARAM( 9) = 0.0D0           ! ORACLE (MIDACO)
            PARAM(10) = 100000D0        ! PARETOMAX (MIDACO)
            PARAM(11) = 0.000001D0      ! EPSILON (MIDACO)
            PARAM(12) = -1.000D0        ! BALANCE => focus only on first objective function (DeltaV, not tt) (MIDACO)
            PARAM(13) = 0.0D0           ! CHARACTER (MIDACO)

            ! Prevent cores racing ahead to optimise over undefined memory

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

        end subroutine variable_init

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Unloads relevant ephemerides; closes currently-used status file.
        !
        ! REVISION HISTORY:
        ! 10 Oct 2019 - Initial version
        ! 03 Jul 2020 - Refactor; Doxygen support
        !
        !--------------------------------------------------------------------------- 

        subroutine intermediate_variable_destruct()

            call UNLOAD(ephemeris_prefix//targ_can//'.bsp')
            close(iunit)

        end subroutine intermediate_variable_destruct

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Initialises all variables needed to reset/restart a MIDACO run
        !
        ! REVISION HISTORY:
        ! 10 Oct 2019 - Initial version
        ! 4 Jul 2020 - Refactor; Doxygen support
        !
        !--------------------------------------------------------------------------- 

        subroutine intermediate_variable_init()

            use midaco_interface
            implicit none

            integer :: temp

            call FURNSH(ephemeris_prefix//targ_can//'.bsp')

            !
            ! Initialise transfer epoch bounds - get the state of the candidate at correct epoch
            !

            call STR2ET('Jan 1 2025 00:00', time_lower)
            call STR2ET('Dec 30 2099 00:00', time_upper)

            XL(1) = time_lower                                                          ! Minimum transfer epoch (ephemeris seconds)
            XL(2) = minimum_transfer_time                                               ! Minimum transfer duration (seconds)
            XL(3) = -25.d0                                                              ! t_end - 1500 days
            XL(4) = 2                                                                   ! n_mnfd
            XL(5) = 2                                                                   ! J
            XU(1) = time_upper                                                          ! Maximum transfer epoch (ephemeris seconds)
            XU(2) = maximum_transfer_time                                               ! Maximum transfer epoch (seconds)
            XU(3) = 0.d0                                                                ! t_end - at the plane
            XU(4) = n_mnfd_disc-1                                                       ! n_mnfd
            XU(5) = num_orbits-1                                                        ! J

            !
            ! Starting point, XOPT
            !

            XOPT = (XL + XU)/2.d0

            !
            ! Maximum function evaluations (MIDACO)
            !

            max_eval = 99999999
            optim_flag = 0
            optim_stop = 0

            ! Printing options

            print_eval = 1000   ! How often to print (MIDACO)
            save_to_file = 0    ! Save solution to text files (MIDACO)

            ! Choose MIDACO parameters (FOR ADVANCED USERS)

            PARAM( 1) = 0.0D0           ! ACCURACY (MIDACO)
            PARAM( 2) = 0.0D0           ! SEED (MIDACO)
            PARAM( 3) = 0.0D0           ! FSTOP (MIDACO)
            PARAM( 4) = 0.0D0           ! ALGOSTOP (MIDACO)
            PARAM( 5) = 0.0D0           ! EVALSTOP (MIDACO)
            PARAM( 6) = 100.D0          ! FOCUS (MIDACO)
            PARAM( 7) = 1000.0D0        ! ANTS (MIDACO)
            PARAM( 8) = 50.0D0          ! KERNEL (MIDACO)
            PARAM( 9) = 0.0D0           ! ORACLE (MIDACO)
            PARAM(10) = 100000D0        ! PARETOMAX (MIDACO)
            PARAM(11) = 0.000001D0      ! EPSILON (MIDACO)
            PARAM(12) = -1.000D0        ! BALANCE => focus only on first objective function (DeltaV, not tt) (MIDACO)
            PARAM(13) = 0.0D0           ! CHARACTER (MIDACO)

            !
            ! Reset work arrays (MIDACO)
            !

            iw = 0
            rw = 0.d0
            pf = 0.d0

            temp = reset_midaco_run()

        end subroutine intermediate_variable_init

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Destroys all constructed MPI variables/pointers/dataspaces used.
        !
        ! REVISION HISTORY:
        ! 10 Oct 2019 - Initial version
        !
        !--------------------------------------------------------------------------- 

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

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Destroys/frees all pointers and data loaded in during the program execution.
        !
        ! REVISION HISTORY:
        ! 20 Mar 2018 - Initial version
        !--------------------------------------------------------------------------- 

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
            call UNLOAD(ephemeris_prefix//targ_can//'.bsp')

            close(iunit)

        end subroutine variable_destruct

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Loads databases into memory.
        !
        ! REVISION HISTORY:
        ! 20 Mar 2018 - Initial version
        ! 04 Jul 2020 - Refactored; Doxygen
        !
        !--------------------------------------------------------------------------- 

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

            if (is_loaded) error STOP "Database has already been loaded. Likely program logic error. Aborting..."

            ! Open file; set number of targets (== number of lines) to zero

            write(*, '(A)', advance="no") " Determining number of orbits in the data-file..."
            open(98, file=original_orbit_data)

            print *, "Number of target points is ", num_targets
            print *, "Number of orbits is ", num_orbits
            print *, "Number of n_mnfd is ", n_mnfd_disc
            print *, "Number of t_end is ",  t_end_disc

            ! Now we know the number of lines, rewind the file pointer and allocate the database
            ! file into memory (likely very large, so approach with caution on low-memory machines)

            rewind(97)
            rewind(98)

            ! Now read the data into memory

            write(*, '(A)', advance="no") " Reading dataset into memory..."

            do k = 1, num_orbits

                do j = 1, n_mnfd_disc

                        read(97, *) dataset(k, j, :)

                       read(98, *) perturbed_conds_dataset(j, k, :)

                end do ! j

            end do ! k

            ! Close file; set is_loaded to true.

            perturbed_conds_dataset = 0.d0

            close(97)
            close(98)

            is_loaded = .true.

        end subroutine load_data

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Load the list of asteroids to optimise into memory.
        !
        ! REVISION HISTORY:
        ! 20 Mar 2018 - Initial version
        ! 04 Jul 2020 - Refactor; Doxygen.
        !
        !--------------------------------------------------------------------------- 

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

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Print a nice message to the screen.
        !
        ! --------------------------------------------------------------------------
        subroutine nice_printing()

            print *, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            print *, " "
            print *, "  WELCOME TO THE FORTRAN ASTEROID RETRIEVAL TOOL (FART)  "
            print *, " "
            print *, "                  Author: Jack Tyler                     "
            print *, "              Email: jack.tyler@soton.ac.uk              "
            print *, " "
            print *, " This version last updated 2020-07-01                "
            print *, " "
            print *, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            print *, " "
            print *, "Parameters specified"
            print *, "~~~~~~~~~~~~~~~~~~~~"
            print *, " "
            print *, "Input data-file: ", datafile
            print *, "Original orbit data:", original_orbit_data
            print *, " "
            print *, "Beginning run "
            print *, "~~~~~~~~~~~~~ "

        end subroutine nice_printing


end module variable_initialisation
        
