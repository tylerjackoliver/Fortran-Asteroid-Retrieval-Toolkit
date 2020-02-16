module variable_initialisation

    use                 precision_kinds
    use                 constants
    use                 state_determination
    use                 problem_parameters
    use                 ancillary_data
    use, intrinsic   :: iso_c_binding, only: c_pointer, c_f_pointer

    integer, parameter                  :: O = 2                                                        ! Number of objectives
    integer, parameter                  :: N = 4                                                        ! Number of variables
    integer, parameter                  :: NI = 0                                                       ! Number of integer variables
    integer, parameter                  :: M = 0                                                        ! Number of contrains
    integer, parameter                  :: ME = 0                                                       ! Number of equality constraints
    integer, parameter                  :: LIW = 5000                                                   ! Integer workspace
    integer, parameter                  :: LRW = 20000                                                  ! Real workspace
    integer, parameter                  :: LPF = 2000000                                                  ! Number of pareto front

    double precision                    :: XL(4), XU(4)                                                 ! Upper and lower bounds
    double precision                    :: XOPT(4)                                                      ! OPtimisation variables
    double precision                    :: F(2)                                                         ! Objectives
    double precision                    :: G(4)                                                         ! Constraint arrays; initialised but never used
    double precision                    :: PARAM(13)                                                    ! MIDACO parameters
    double precision                    :: RW(LRW), PF(LPF)                                             ! Workspace and pareto front

    double precision, allocatable       :: pareto_front_data(:, :)                                ! Stores the custom pareto_front_data for generation of ancillary data

    integer                             :: optim_flag                                                   ! Optimiser information flag
    integer                             :: optim_stop                                                   ! Optimiser stopping variable
    integer                             :: IW(LIW)                                                      ! Integer workspace
    integer                             :: max_time                                                     ! Maximum walltime for optimisations
    integer                             :: max_eval                                                     ! Maximum function evaluations
    integer                             :: print_eval                                                   ! How often to print
    integer                             :: save_to_file                                                 ! Output verbosity
    integer                             :: number_of_pareto_solutions                                   ! ...number of Pareto solutions

    character*60                        :: key                                                          ! License key bit

    double precision                    :: time_lower                                                   ! Lower bound for the optimisation time - dummy here
    double precision                    :: time_upper                                                   ! Upper bound for the optimisation time - dummy here

    logical                             :: is_loaded = .false.

    integer                             :: num_targets                                                  ! Number of items in the input file
    integer                             :: num_orbits                                                   ! Number of orbits in the file

    integer, parameter                  :: iunit = 70                                                       ! Unit for the Pareto file
    
    integer, parameter                  :: t_end_disc = 100                                             ! Parameterisation in backwards time
    integer, parameter                  :: n_mnfd_disc = 360                                            ! Parameterisation around the orbit

    ! MPI variables

    integer                             :: mpi_err                                                      ! MPI error integer
    integer                             :: mpi_id_world                                                 ! Core ID of the MPI process
    integer                             :: mpi_id_shared                                                ! ID of the process in the shared-space
    integer                             :: mpi_world_size 
    integer                             :: mpi_name_string_size
    integer                             :: mpi_color

    ! MPI communicators

    integer(kind=MPI_COMM)              :: node_communicator                                            ! Intra-node communicator
    integer(kind=MPI_COMM)              :: node_master_communicator                                     ! Master nodes on each core

    ! MPI node naming

    character(MPI_MAX_PROCESSOR_NAME)   :: proc_name

    ! MPI window variables

    TYPE(C_PTR)                         :: dataset_pointer                                              ! C-style pointer to the shared memory array
    TYPE(C_PTR)                         :: perturbed_conds_pointer                                      ! C-style pointer to the perturbed conds array
    TYPE(MPI_Win)                       :: dataset_window                                               ! Window that MPI opens to the shared memory
    TYPE(MPI_Win)                       :: perturbed_conds_window                                       ! Window that MPI opens to the shared memory

    integer(kind=MPI_ADDRESS_KIND)      :: dataset_size_bytes
    integer(kind=MPI_ADDRESS_KIND)      :: perturbed_conds_size_bytes
    real*8                              :: dataset_double
    real*8                              :: perturbed_conds_double

    contains

        subroutine mpi_variable_init()

            ! Initialise MPI

            call MPI_INIT(mpi_err)
            call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_id_world, mpi_err)
            call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_world_size, mpi_err)
            call MPI_GET_PROCESSOR_NAME(proc_name, mpi_name_string_size, mpi_err)

#ifdef DEBUG

            print *, "MPI process ID", mpi_id_world, "is running"

#endif

            if (mpi_id_world .eq. 0) then

                print *, "Importing data on rank 0..."
                call load_data()
                print *, "done."

            end if

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

            ! First, broadcast the number of orbits

            call MPI_BCAST(num_orbits, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

            ! Split the global communicator into smaller communicators,
            ! each of which defines a shared-memory region (each node)

            call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, node_communicator, mpi_err)

#ifdef DEBUG

            print *, "I am global rank", mpi_id_world, "rank", mpi_id_shared, "on proc", proc_name

#endif

            ! Now create another communicator for *only* the procs that are masters
            ! on their relative nodes

            if (mpi_id_shared .eq. 0) then

                mpi_color = 1 ! Want this to be in the communicator

            else

                mpi_color = MPI_UNDEFINED ! Do not want this to be in the communicator

            end if

            call MPI_COMM_SPLIT(MPI_COMM_WORLD, mpi_color, 0, node_master_communicator, mpi_err)

            !
            ! Now we have three communicators: WORLD, intra-node, and node-masters
            ! We need to have every node-master have a copy of the target data, so broadcast
            ! the target array *solely* to the node-masters
            !

            ! Allocate memory only for the node masters but NOT the one that has already allocated data

            dataset_size_bytes = num_orbits * t_end_disc * n_mnfd_disc * sizeof(dataset_double)
            perturbed_conds_size_bytes = num_orbits * t_end_disc * n_mnfd_disc * sizeof(perturbed_conds_double)

            ! Allocate memory for access by all on the same node (node_communicator)

            ! Do dataset first...

            call MPI_WIN_ALLOCATE_SHARED(dataset_size_bytes, 1, MPI_INFO_NULL, node_communicator, dataset_pointer, dataset_window, mpi_err)

            ! ...and now the perturbed_conds array

            call MPI_WIN_ALLOCATE_SHARED(perturbed_conds_size_bytes, 1, MPI_INFO_NULL, node_communicator, perturbed_conds_pointer, perturbed_conds_window, mpi_err)

            ! Obtain the locaiton of the shared memory segments for non-master cores
            !
            ! Call the same name so that the syntax is the same
            !

            if (mpi_id_shared .ne. 0) then

                call MPI_WIN_SHARED_QUERY(dataset_window, mpi_id_shared, dataset_size_bytes, 1, dataset_pointer, mpi_err)
                call MPI_WIN_SHARED_QUERY(perturbed_conds_window, mpi_id_shared, perturbed_conds_size_bytes, 1, perturbed_conds_pointer, mpi_err)

            end if

            ! Now convert the pointer to an array

            call C_F_POINTER(dataset_pointer, dataset, (/num_orbits, t_end_disc, n_mnfd_disc, 7/))
            call C_F_POINTER(perturbed_conds_pointer, perturbed_conds, (/num_orbits, n_mnfd_disc, 7/))

            ! Initialise MPI RMA calls

            call MPI_WIN_START(dataset_window, mpi_err)
            call MPI_WIN_START(perturbed_conds_window, mpi_err)

        end subroutine mpi_variable_init

        subroutine variable_init()

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
            XL(2) = 1.D0 * 86400.D0                                                        ! Minimum transfer duration (seconds)
            XL(3) = 1                                                                      ! t_end
            XL(4) = 1                                                                      ! n_mnfd
            XU(1) = time_upper                                                             ! Maximum transfer epoch (ephemeris seconds)
            XU(2) = 1500.D0 * 86400.D0                                                     ! Maximum transfer epoch (seconds)
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

        subroutine intermediate_variable_destruct()

            call UNLOAD('../data/'//targ_can//'.bsp')

        end subroutine intermediate_variable_destruct

        subroutine intermediate_variable_init()

            call FURNSH('../data/'//targ_can//'.bsp')

            ! Initialise transfer epoch bounds - get the state of the candidate at correct epoch

            call GET_STATE(targ_can, time_lower, time_upper)

            ! Optimiser bounds

            XL(1) = time_lower                                                             ! Minimum transfer epoch (ephemeris seconds)
            XL(2) = 1.D0 * 86400.D0                                                        ! Minimum transfer duration (seconds)
            XL(3) = 1                                                                      ! t_end
            XL(4) = 1                                                                      ! n_mnfd
            XU(1) = time_upper                                                             ! Maximum transfer epoch (ephemeris seconds)
            XU(2) = 1500.D0 * 86400.D0                                                     ! Maximum transfer epoch (seconds)
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

        end subroutine intermediate_variable_init

        subroutine mpi_variable_destruct()

            ! ////////////////////////////////////////////////////////////
            !
            ! Ensures graceful program closure
            !
            ! ////////////////////////////////////////////////////////////

            deallocate(dataset)
            deallocate(perturbed_conds)

            ! Unload SPICE kernels
                    
            call UNLOAD('../data/de414.bsp')                    
            call UNLOAD('../data/naif0008.tls')
            call UNLOAD('../data/'//targ_can//'.bsp')

            ! Fence-off access to the RMA space

            call MPI_WIN_FENCE(0, dataset_window, mpi_err)
            call MPI_WIN_FENCE(0, perturbed_conds_window, mpi_err)

            ! Finalise access

            call MPI_WIN_COMPLETE(dataset_window, mpi_err)
            call MPI_WIN_COMPLETE(perturbed_conds_window, mpi_err)

            ! Free the space

            call MPI_WIN_FREE(dataset_window, mpi_err)
            call MPI_WIN_FREE(perturbed_conds_window, mpi_err)

            call MPI_FINALIZE(mpi_err)

            close(iunit)

        end subroutine mpi_variable_destruct

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

            write(*, '(A)', advance="no") "Determining number of lines..."
            num_targets = 0
            open(97, file=datafile)
            open(98, file="../data/2020-01-28_L2PlanarPertubedConds.csv")

            ! Read through in do-loop until EOF to determine number of lines

            do

                read(97,*,iostat=io_state)                                                  ! Not interested in contents yet
                
                if (io_state .lt. 0) then                                                   ! If an error occured (likely EOF)

                    exit                                                                    ! Do loop

                end if

                num_targets = num_targets + 1

            end do

            write(*, '(A)') "done."

            ! Now we know the number of lines, rewind the file pointer and allocate the database
            ! file into memory (likely very large, so approach with caution on low-memory machines)

            rewind(97)

            ! Move from number of targets -> number of orbits

            num_targets = num_targets

            num_orbits = num_targets / (n_mnfd_disc * t_end_disc)

            ! Now read the data into memory

            write(*, '(A)', advance="no") "Reading dataset into memory..."

            do k = 1, num_orbits

                do j = 1, n_mnfd_disc

                    do i = 1, t_end_disc

                        read(97, *) dum, dataset(i, j, k, :) ! Dum == time

                    end do ! i

                    read(98) perturbed_conds(j, k, :) ! Dum == time

                end do ! j

            end do ! k

            write(*, '(A, I8)') "done. Number of targets is", num_targets
            print '(A, I8)', "Number of orbits is", num_orbits
            print '(A, I8)', "Number of n_mnfd is", n_mnfd_disc
            print '(A, I8)', "Number of t_end is",  t_end_disc

            ! Close file; set is_loaded to true.

            close(97)
            close(98)

            is_loaded = .true.

        end subroutine load_data

        subroutine get_pareto_front() 

            double precision    :: current_working_state(4)
            double precision    :: final_velocity

            integer             :: orbit_choice
            integer             :: pareto_solution
            integer             :: variable_iter
            integer             :: orbit_num

            double precision    :: r1(6)
            double precision    :: v1(3)
            double precision    :: r2(6)
            double precision    :: v2(3)
            double precision    :: target_point(6)
            double precision    :: tmani

            open(iunit, file="../data/paretoFront_"//targ_can)

            ! Get number of Pareto solutions

            number_of_pareto_solutions = PF(1) ! First element stores the number of Pareto solutions

            ! Now we know the number of Pareto solution, we can get the solutions

            call write_pareto_head()

            do pareto_solution = 1, number_of_pareto_solutions

                final_velocity = PF(2 + O*(pareto_solution-1))

                do variable_iter = 1, N

                    current_working_state(variable_iter) = PF(2 + O * PARETOMAX + M * PARETOMAX + &
                                                              n*(pareto_solution-1) + variable_iter - 1)

                end do

                call generate_ancillary_data(current_working_state(1), current_working_state(2), current_working_state(3), &
                                             current_working_state(4), orbit_num, r1, v1, r2, v2, tmani, target_point)

                call write_pareto_body(final_velocity, current_working_state, orbit_num, r1, v1, r2, v2, target_point, tmani)

            end do

            close(iunit)

        end subroutine get_pareto_front

        subroutine write_pareto_head()

            write(iunit, '(A)') "## Pareto front solution file "
            write(iunit, '(A)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            write(iunit, '(A)') " "
            write(iunit, '(A)') "deltaV \t t0 \t tt \t tend \t nmnfd \t J \t v1x v1y v1z [km/s] \t v2x v2y v2z [km/s] \t tmani [rS] \t xtarg [CR3BP]"

        end subroutine write_pareto_head

        subroutine write_pareto_body(final_velocity, current_working_state, orbit_num, v1, v2, target_point, tmani)

            double precision, intent(in)    :: current_working_state(4)
            double precision, intent(in)    :: final_velocity
            double precision, intent(in)    :: tmani

            integer, intent(in)             :: orbit_choice
            integer, intent(in)             :: pareto_solution
            integer, intent(in)             :: variable_iter
            integer, intent(in)             :: orbit_num

            double precision, intent(in)    :: v1(3)
            double precision, intent(in)    :: v2(3)
            double precision, intent(in)    :: target_point(6)
            
            write(iunit, *) final_velocity, current_working_state, orbit_choice, v1, v2, tmani, target_point

        end subroutine write_pareto_body

end module variable_initialisation
        
