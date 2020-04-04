module variable_initialisation

    use                 precision_kinds
    use                 mpi
    use                 constants
    use                 state_determination
    use                 problem_parameters
    use, intrinsic   :: iso_c_binding

    integer, parameter                  :: O = 2                                                        ! Number of objectives
    integer, parameter                  :: N = 5                                                        ! Number of variables
    integer, parameter                  :: NI = 0                                                       ! Number of integer variables
    integer, parameter                  :: M = 0                                                        ! Number of contrains
    integer, parameter                  :: ME = 0                                                       ! Number of equality constraints
    integer, parameter                  :: LIW = 5000                                                   ! Integer workspace
    integer, parameter                  :: LRW = 20000                                                  ! Real workspace
    integer, parameter                  :: LPF = 2000000                                                  ! Number of pareto front

    double precision                    :: XL(5), XU(5)                                                 ! Upper and lower bounds
    double precision                    :: XOPT(5)                                                      ! OPtimisation variables
    double precision                    :: F(2)                                                         ! Objectives
    double precision                    :: G(4)                                                         ! Constraint arrays; initialised but never used
    double precision                    :: PARAM(13)                                                    ! MIDACO parameters
    double precision                    :: RW(LRW), PF(LPF)                                             ! Workspace and pareto front

    double precision, allocatable       :: pareto_front_data(:, :)                                      ! Stores the custom pareto_front_data for generation of ancillary data

    integer                             :: optim_flag                                                   ! Optimiser information flag
    integer                             :: optim_stop                                                   ! Optimiser stopping variable
    integer                             :: IW(LIW)                                                      ! Integer workspace
    integer                             :: max_time                                                     ! Maximum walltime for optimisations
    integer                             :: max_eval                                                     ! Maximum function evaluations
    integer                             :: print_eval                                                   ! How often to print
    integer                             :: save_to_file                                                 ! Output verbosity
    integer                             :: number_of_pareto_solutions                                   ! ...number of Pareto solutions
    integer                             :: iostate
    integer                             :: num_candidates
    integer                             :: candidates_per_processor
    integer                             :: host

    integer, allocatable                :: candidates_to_consider(:)

    character*60                        :: key                                                          ! License key bit

    double precision                    :: time_lower                                                   ! Lower bound for the optimisation time - dummy here
    double precision                    :: time_upper                                                   ! Upper bound for the optimisation time - dummy here

    logical                             :: is_loaded = .false.

    integer                             :: num_targets                                                  ! Number of items in the input file
    integer                             :: num_orbits                                                   ! Number of orbits in the file
    integer(8)                          :: num_orbits_dble

    Integer                             :: iunit                                                        ! Unit for the Pareto file

    ! MPI variables

    integer                             :: mpi_err                                                      ! MPI error integer
    integer                             :: mpi_id_world                                                 ! Core ID of the MPI process
    integer                             :: mpi_id_shared                                                ! ID of the process in the shared-space
    integer                             :: mpi_world_size 
    integer                             :: mpi_name_string_size
    integer                             :: mpi_color

    ! MPI communicators

    integer                             :: node_communicator                                            ! Intra-node communicator
    integer                             :: node_master_communicator                                     ! Master nodes on each core

    ! MPI node naming

    character(MPI_MAX_PROCESSOR_NAME)   :: proc_name

    ! MPI window variables

    TYPE(C_PTR)                         :: dataset_pointer                                              ! C-style pointer to the shared memory array
    TYPE(C_PTR)                         :: perturbed_conds_pointer                                      ! C-style pointer to the perturbed conds array
    integer                             :: dataset_window                                               ! Window that MPI opens to the shared memory
    integer                             :: perturbed_conds_window                                       ! Window that MPI opens to the shared memory

    integer                             :: info

    integer(kind=MPI_ADDRESS_KIND)      :: dataset_size_bytes
    integer(kind=MPI_ADDRESS_KIND)      :: perturbed_conds_size_bytes
    integer(8)                          :: row_length
    integer(8)                          :: row_length_perturbed
    integer(8)                          :: double_size
    integer                             :: disp_unit
    integer(8)                          :: t_end_disc_dble = t_end_disc
    integer(8)                          :: n_mnfd_disc_dble = n_mnfd_disc
    integer(8)                          :: to_send_size
    real*8                              :: dataset_double
    real*8                              :: perturbed_conds_double

    contains

        subroutine mpi_variable_init()

            integer                     :: i

            disp_unit=1

            ! Initialise MPI

            call MPI_INIT(mpi_err)
            call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_id_world, mpi_err)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_world_size, mpi_err)
            call MPI_GET_PROCESSOR_NAME(proc_name, mpi_name_string_size, mpi_err)

            if (mpi_id_world .eq. 0) then

                call nice_printing()
                call determine_num_orbits()

            end if

            ! First, broadcast the number of orbits

#ifdef DEBUG

            if (mpi_id_world .eq. 0) then 

                print *, "Sending number of orbits, ", num_orbits

            end if

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

#endif

            call MPI_BCAST(num_orbits, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

#ifdef DEBUG

            print *, "I am rank", mpi_id_world, "and I have the number of orbits as", num_orbits

#endif

            num_orbits_dble = num_orbits

            ! Split the global communicator into smaller communicators,
            ! each of which defines a shared-memory region (each node)

            call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_communicator, mpi_err)
            call MPI_COMM_RANK(node_communicator, mpi_id_shared, mpi_err)

#ifdef DEBUG

            print *, "I am global rank", mpi_id_world, "rank", mpi_id_shared, "on proc", proc_name

#endif

            ! Now create another communicator for *only* the procs that are masters
            ! on their relative nodes

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

            row_length = 7
            row_length_perturbed = 6
            double_size = sizeof(dataset_double)
            
            if (mpi_id_shared .eq. 0) then

                dataset_size_bytes = int(num_orbits * t_end_disc * n_mnfd_disc * row_length * double_size, MPI_ADDRESS_KIND)
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

            ! Now convert the C pointer to a Fortran pointer (syntax for accessing elements is the same as standard ALLOCATEd arrays)
             
            call C_F_POINTER(dataset_pointer, dataset, (/t_end_disc_dble, n_mnfd_disc_dble, &
                                                        num_orbits_dble, row_length/))
            call C_F_POINTER(perturbed_conds_pointer, perturbed_conds_dataset, &
                (/n_mnfd_disc_dble, num_orbits_dble, row_length_perturbed/))
  
            ! To avoid racing on the IO, have the global master _only_ read in the data

            if (mpi_id_world .eq. 0) then

                call load_data()
                call load_targets()

            end if

            ! Tell all nodes to WAIT PLEASE, BACK BEHIND THE LINE

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

            !
            ! Now the global master should BCAST all of the data it just read-in to other node-masters
            ! However, the data it is sending overflows the internal C-type `int` used in MPI, so we can only
            ! send up to 2**32-1 bits at a time.
            !
            ! Here, we eat the overheads and send 'stripes' of the dataset and perturbed_conds dataset out to
            ! the other node-masters. This will not overflow *our* dataset, and makes reassembling the data on the other end
            ! far easier compared to using e.g. MPI structures (and to be honest, I cba to implement MPI STRUCTURES today)
            !
            ! Provided all workers call this subroutine, then it is also blocking/synchronous
            ! 

            if (mpi_id_shared .eq. 0) then

#ifdef DEBUG

                print *, "Rank ", mpi_id_world, "is sending data"

#endif

                host = 0

                do i = 1, t_end_disc

                    ! Send dataset

                    to_send_size = size(dataset(i, :, :,:))
                    call MPI_BCAST(dataset(i, :, :, :), to_send_size, MPI_DOUBLE_PRECISION, host, node_master_communicator, mpi_err)
                    
                end do
                
                print *, "Rank", mpi_id_world, n_mnfd_disc

                do i = 1, n_mnfd_disc

                        ! Send perturbed_conds

                    to_send_size = size(perturbed_conds_dataset(i, :, :))
                    call MPI_BCAST(perturbed_conds_dataset(i,:,:), to_send_size, MPI_DOUBLE_PRECISION, host,&
                                         node_master_communicator, mpi_err) 

                end do

            end if
                
            ! Send perturbed_conds_dataset


#ifdef DEBUG

            print *, "rank", mpi_id_world, "has finished sending/receiving data"
            
#endif

            !
            ! Tell all nodes to wait for the node masters to finish writing their love letters to each other,
            ! or break their hearts in a vice after running off with the dog ISN'T THAT RIGHT SCARLETT
            !

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

            !
            ! Now send the number of candidates that we have to all of the other nodes
            !

            host = 0
            call MPI_BCAST(num_candidates, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)    
            
            if (mpi_id_world .ne. 0) allocate(targ_can_array(num_candidates))

            !
            ! And use that information to allocate arrays with the desired targets in
            !

            !
            ! ...and populate them with the candidates they should be looking at
            !

            to_send_size = size(targ_can_array)
            call MPI_BCAST(targ_can_array, to_send_size, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
#ifdef DEBUG

            print *, "Rank: ", mpi_id_world, "has that bit o the dataset as", dataset(864, 90, 2, :), &
                " and the bits before and after ", dataset(864, 90, 1, :), &
                dataset(864, 90, 3, :)

#endif

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

!            num_orbits = num_targets / (n_mnfd_disc * t_end_disc)
            num_orbits = 500
            num_targets = num_orbits * n_mnfd_disc * t_end_disc

        end subroutine determine_num_orbits

        subroutine variable_init()


            ! Load SPICE kernels

            iunit = 151+mpi_id_world
            call FURNSH('../data/de414.bsp')                    
            call FURNSH('../data/naif0008.tls')
            call FURNSH('../data/'//targ_can//'.bsp')
            open(iunit, file="../data/paretoFront_"//targ_can)

            ! Initialise transfer epoch bounds - get the state of the candidate at correct epoch

            call GET_STATE(targ_can, time_lower, time_upper)

            ! Optimiser bounds

            XL(1) = time_lower                                                             ! Minimum transfer epoch (ephemeris seconds)
            XL(2) = 1.D0 * 86400.D0                                                        ! Minimum transfer duration (seconds)
            XL(3) = maxval(dataset(1000,:,:,1))                                            ! t_end - should be the *least negative* backwards integration time of all orbits
            XL(4) = 1                                                                      ! n_mnfd
            XL(5) = 1                                                                      ! J
            XU(1) = time_upper                                                             ! Maximum transfer epoch (ephemeris seconds)
            XU(2) = 1500.D0 * 86400.D0                                                     ! Maximum transfer epoch (seconds)
            XU(3) = maxval(dataset(1,:,:,1))                                               ! t_end - should be the *most negative* time of all pi/8 planes
            XU(4) = n_mnfd_disc                                                            ! n_mnfd
            XU(5) = num_orbits                                                                    ! J


            ! Starting point, XOPT

            XOPT = (/(time_lower+time_upper)*.5d0, 750.d0 * 86400d0, (XL(3)+XU(3))/2.d0, &
                    n_mnfd_disc/2.d0, floor(num_orbits/2.d0) * 1.d0/) ! EXACT MIDDLE OF THE SET

            print *, "Initial point", XOPT
            
            ! Maximum function evaluations

            max_eval = 99999999
            max_time = 300 ! 60

            ! Printing options

            print_eval = 100
            save_to_file = 0 ! Save solution to text files

            ! Choose MIDACO parameters (FOR ADVANCED USERS)

            PARAM( 1) = 0.0D0           ! ACCURACY
            PARAM( 2) = 0.0D0           ! SEED
            PARAM( 3) = 0.0D0           ! FSTOP
            PARAM( 4) = 0.0D0           ! ALGOSTOP
            PARAM( 5) = 0.0D0           ! EVALSTOP
            PARAM( 6) = 0.D0            ! FOCUS
            PARAM( 7) = 0.0D0           ! ANTS
            PARAM( 8) = 0.0D0           ! KERNEL
            PARAM( 9) = 0.0D0           ! ORACLE
            PARAM(10) = 100000D0        ! PARETOMAX
            PARAM(11) = 0.000001D0      ! EPSILON  
            PARAM(12) = -1.000D0        ! BALANCE => focus only on first objective function (DeltaV, not tt)
            PARAM(13) = 0.0D0           ! CHARACTER 

            ! If the dataset has not been loaded, AND the core is a master core on any node (rank zero on its own
            ! communicator), then load the dataset in here using the shared memory space as a target

            ! if (.not. is_loaded .and. mpi_id_shared .eq. 0) call load_data()

            ! Prevent cores racing ahead to optimise over undefined memory

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

        end subroutine variable_init

        subroutine intermediate_variable_destruct()

            call UNLOAD('../data/'//targ_can//'.bsp')

        end subroutine intermediate_variable_destruct

        subroutine intermediate_variable_init()

            use midaco_interface
            implicit none

            integer :: temp

            call FURNSH('../data/'//targ_can//'.bsp')

            ! Initialise transfer epoch bounds - get the state of the candidate at correct epoch

            call GET_STATE(targ_can, time_lower, time_upper)
            
            ! Optimiser bounds

            XL(1) = time_lower                                                             ! Minimum transfer epoch (ephemeris seconds)
            XL(2) = 1.D0 * 86400.D0                                                        ! Minimum transfer duration (seconds)
            XL(3) = maxval(dataset(1000,:,:,1))                                            ! t_end - should be the *least negative* backwards integration time of all orbits
            XL(4) = 1                                                                      ! n_mnfd
            XL(5) = 1                                                                      ! J
            XU(1) = time_upper                                                             ! Maximum transfer epoch (ephemeris seconds)
            XU(2) = 1500.D0 * 86400.D0                                                     ! Maximum transfer epoch (seconds)
            XU(3) = maxval(dataset(1,:,:,1))                                               ! t_end - should be the *most negative* time of all pi/8 planes
            XU(4) = n_mnfd_disc                                                            ! n_mnfd
            XU(5) = num_orbits                                                                    ! J

            ! Starting point, XOPT

            XOPT = (/(time_lower+time_upper)*.5d0, 750.d0 * 86400d0, (XL(3)+XU(4)), &
                    n_mnfd_disc/2.d0, floor(num_orbits/2.d0) * 1.d0/) ! EXACT MIDDLE OF THE SET

            ! Maximum function evaluations

            max_eval = 99999999
            max_time = 10 ! 2 days
            optim_flag = 0
            optim_stop = 0

            ! Printing options

            print_eval = 10000
            save_to_file = 0 ! Save solution to text files

            ! Choose MIDACO parameters (FOR ADVANCED USERS)

            PARAM( 1) = 0.0D0           ! ACCURACY
            PARAM( 2) = 0.0D0           ! SEED
            PARAM( 3) = 0.0D0           ! FSTOP
            PARAM( 4) = 0.0D0           ! ALGOSTOP
            PARAM( 5) = 0.0D0           ! EVALSTOP
            PARAM( 6) = 0.D0            ! FOCUS
            PARAM( 7) = 0.0D0           ! ANTS
            PARAM( 8) = 0.0D0           ! KERNEL
            PARAM( 9) = 0.0D0           ! ORACLE
            PARAM(10) = 100000D0        ! PARETOMAX
            PARAM(11) = 0.000001D0      ! EPSILON  
            PARAM(12) = -1.000D0        ! BALANCE => focus only on first objective function (DeltaV, not tt)
            PARAM(13) = 0.0D0           ! CHARACTER 

            iw = 0
            rw = 0.d0
            pf = 0.d0

            temp = reset_midaco_run()

        end subroutine intermediate_variable_init

        subroutine mpi_variable_destruct()

            ! ////////////////////////////////////////////////////////////
            !
            ! Ensures graceful program closure
            !
            ! ////////////////////////////////////////////////////////////

            ! Unload SPICE kernels
                    
            call UNLOAD('../data/de414.bsp')                    
            call UNLOAD('../data/naif0008.tls')
            call UNLOAD('../data/'//targ_can//'.bsp')

            ! Free the shared-memory space

            call MPI_WIN_FREE(dataset_window, mpi_err)
            call MPI_WIN_FREE(perturbed_conds_window, mpi_err)

            call MPI_FINALIZE(mpi_err)

            close(iunit)

        end subroutine mpi_variable_destruct

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

            if (is_loaded) error STOP "Database has already been loaded. Likely program logic error. Aborting..."

            ! Open file; set number of targets (== number of lines) to zero

            write(*, '(A)', advance="no") " Determining number of orbits in the data-file..."
            num_targets = 0
            open(98, file=original_orbit_data)

            ! ! Move from number of targets -> number of orbits

            ! num_orbits = num_targets / (n_mnfd_disc * t_end_disc)

            ! write(*, *) "done."
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

                    do i = 1, t_end_disc

                        read(97, *) dataset(i, j, k, :) ! Dum == time

                    end do ! i
                        
                    read(98, *) perturbed_conds_dataset(j, k, :) ! Dum == time

                end do ! j

            end do ! k

            ! Close file; set is_loaded to true.

            close(97)
            close(98)

            is_loaded = .true.

        end subroutine load_data


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
            print *, " This version last updated 2020-02-03                    "
            print *, " "
            print *, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            print *, " "
            print *, "Parameters specified"
            print *, "~~~~~~~~~~~~~~~~~~~~"
            print *, " "
            print *, "Input data-file: ", datafile
            print *, "Original orbit data:", original_orbit_data
            print *, "Target candidate: ", targ_can
            print *, " "
            print *, "Beginning run "
            print *, "~~~~~~~~~~~~~ "

        end subroutine nice_printing


end module variable_initialisation
        
