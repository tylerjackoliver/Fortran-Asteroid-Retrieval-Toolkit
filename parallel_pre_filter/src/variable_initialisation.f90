module variable_init

    use precision_kinds
    use constants
    use iso_c_binding, only : c_ptr, c_f_pointer
    use problem_parameters
    use utility_functions
    use mpi

    implicit none
    
    integer, allocatable                        :: database(:)                                                  ! SPK ID Database
    integer                                     :: candidate_count = 0                                          ! Number of candidates found
    integer                                     :: itercount = 0                                                ! Main loop sentinel
    integer                                     :: num_orbits = 0
    integer                                     :: num_weeks = 0
    integer                                     :: num_candidates = 0

    double precision                            :: dum                                                          ! Dummy variable

    ! MPI variables

    integer                                     :: mpi_err                                                      ! MPI error integer
    integer                                     :: mpi_id_world                                                 ! Core ID of the MPI process
    integer                                     :: mpi_id_shared                                                ! ID of the process in the shared-space
    integer                                     :: mpi_world_size 
    integer                                     :: mpi_name_string_size
    integer                                     :: mpi_color
    integer                                     :: double_size
    integer                                     :: host
    integer                                     :: row_length

    ! MPI communicators

    integer                                     :: node_communicator                                            ! Intra-node communicator
    integer                                     :: node_master_communicator                                     ! Master nodes on each core

    ! MPI node naming

    character(MPI_MAX_PROCESSOR_NAME)           :: proc_name

    ! MPI window variables

    TYPE(C_PTR)                                 :: manifold_data_elements_ptr                                   ! C-style pointer to the shared memory array
    integer                                     :: manifold_elements_window                                     ! Window that MPI opens to the shared memory

    integer                                     :: info

    integer(kind=MPI_ADDRESS_KIND)              :: dataset_size_bytes
    integer                                     :: disp_unit
    integer                                     :: to_send_size

    double precision                            :: dataset_double = 1.d0                                        ! To measuresize of double precision for MPI

    contains

        subroutine mpi_variable_initialisation()

            use precision_kinds
            use constants

            implicit none

            integer                             :: iostate
            integer                             :: i

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
#ifdef VERBOSE

            if (mpi_id_world .eq. 0) then 

                print *, "Sending number of orbits, ", num_orbits

            end if
#endif
            call MPI_BCAST(num_orbits, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

            ! Split the global communicator into smaller communicators,
            ! each of which defines a shared-memory region (each node)

            call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_communicator, mpi_err)
            call MPI_COMM_RANK(node_communicator, mpi_id_shared, mpi_err)

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

            row_length = 8
            double_size = sizeof(dataset_double)
            
            if (mpi_id_shared .eq. 0) then

                dataset_size_bytes = int(num_orbits * row_length * double_size, MPI_ADDRESS_KIND)

            else

                dataset_size_bytes = 0

            end if

            call MPI_WIN_ALLOCATE_SHARED(dataset_size_bytes, disp_unit, MPI_INFO_NULL, node_communicator, &
                                         manifold_data_elements_ptr, manifold_elements_window, mpi_err)

            !
            ! If we're not the master thread, then we won't know the location of the memory space; thus, query the location of the memory here for use later
            !

            if (mpi_id_shared /= 0) then

                call MPI_WIN_SHARED_QUERY(manifold_elements_window, 0, dataset_size_bytes, disp_unit, &
                                        manifold_data_elements_ptr, mpi_err)

            end if

            ! Now convert the C pointer to a Fortran pointer (syntax for accessing elements is the same as standard ALLOCATEd arrays)
             
            call C_F_POINTER(manifold_data_elements_ptr, manifold_data_elements, (/num_orbits, row_length/))

            ! To avoid racing on the IO, have the global master _only_ read in the data

            if (mpi_id_world .eq. 0) then

                call load_dataset()

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

                host = 0

                do i = 1, num_orbits

                    ! Send dataset

                    to_send_size = size(manifold_data_elements(i,:))
                    call MPI_BCAST(manifold_data_elements(i, :), to_send_size, MPI_DOUBLE_PRECISION, host, &
                                   node_master_communicator, mpi_err)
                    
                end do

            end if

            ! Tell all nodes to WAIT PLEASE, BACK BEHIND THE LINE

            call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

            !
            ! Now load the week IDs to consider into memory
            !

            if (mpi_id_world .eq. 0) then

                call load_week_ids()

            end if

            !
            ! Broadcast the week IDs to all other workers
            !

            ! First, send num_weeks

            call MPI_BCAST(num_weeks, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

            if (mpi_id_world .ne. 0) then

                allocate(dataset(num_weeks))

            end if

            to_send_size = size(dataset)
            call MPI_BCAST(dataset, to_send_size, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

            !
            ! Now broadcast the number of candidates and the candidate database
            !

            if (mpi_id_world .eq. 0) then

                call load_asteroid_database()

            end if

            call MPI_BCAST(num_candidates, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

            if (mpi_id_world .ne. 0) then

                allocate(asteroid_database(num_candidates))

            end if

            to_send_size = size(asteroid_database)
            call MPI_BCAST(asteroid_database, to_send_size, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

            ! Load ephemerides

            call FURNSH(planetary_ephemeris)
            call FURNSH(leapseconds_file)

        end subroutine mpi_variable_initialisation

        subroutine load_dataset()

            integer         :: i

            real(kind=dp)   :: temp(7)

            do i = 1, num_orbits

                read(99, *) temp
                ! print *, "temp", temp
            
                !
                ! Now we have the _synodic_, non-dimensional state, we need to convert to the global inertial frame
                ! in order to determine the correct orbital elements
                !
                ! Note that we DON'T need to globalise since we can assume that mu=1 gives the correct (non-dim) a in the
                ! global frame - just need to multiply by au once we get into the main function 
                !
                ! The epoch for global_rotate can be assumed to be zero, since the magnitude doesn't change (a, e, i) when
                ! we select a different time, and these are the only OEs we're using here
                !

                call GLOBAL_ROTATE(temp(2:7), 0.0d0, manifold_data_state)
            
                !
                ! Epoch may be zero here as we're only interested in orbital element magnitudes (a, e, i) -- not the phasing!
                !
            
                call OSCELT(manifold_data_state, 0.d0, 1.d0, manifold_data_elements(i,:))
            
            end do
            
            close(99)

        end subroutine load_dataset

        subroutine load_week_ids()

            integer :: iostate
            integer :: i

            ! Determine the number of asteroids we're targeting

            write(*, '(A)', advance='no') " Determining number of lines in the SPK dataset file..."
        
            open(unit=98, file=week_database)                                             ! Database file
            
            do
        
                read(98, *, iostat=iostate)
        
                if (iostate .ne. 0) then
        
                    exit
        
                else
        
                    num_weeks = num_weeks + 1

                end if
        
            end do
        
            rewind(98)

            ! Allocate correct size of database

            allocate(dataset(num_weeks))

            write(*,'(A)', advance='no') "done. Reading dataset into memory..."
        
            do i = 1, num_weeks
            
                read(98, *) dataset(i)
            
            end do
        
            close(98)

            print '(A)', " done."

        end subroutine load_week_ids

        subroutine variable_destruct()

            close(output_unit)
            call UNLOAD(planetary_ephemeris)
            call UNLOAD(leapseconds_file)

        end subroutine variable_destruct

        subroutine determine_num_orbits()

            integer :: iostate

        ! Determine number of lines in the target dataset file

            write(*, '(A)', advance='no') " Determining number of lines in the target dataset file..."

            open(unit=99, file=target_dataset)                                                 ! Target OEs
        
            do
        
                read(99, *, iostat=iostate)
        
                if (iostate .ne. 0) then
        
                    exit
        
                else
        
                    num_orbits = num_orbits + 1
        
                end if

            end do
        
            rewind(99)

            print '(A)', "done."

        end subroutine determine_num_orbits

        !
        ! Load in the SPK IDs of the database
        !

        subroutine load_asteroid_database()

            integer :: i
            integer :: iostate

            open(unit=69, file=spk_database)

            do

                read(69, *, iostat=iostate)

                if (iostate .ne. 0) then

                    exit

                end if

                num_candidates = num_candidates + 1

            end do

            ! Now we know the number of lines, allocate and read into memory

            rewind(69)

            allocate(asteroid_database(num_candidates))

            do i = 1, num_candidates

                read(69, *), asteroid_database(i)

            end do

            close(69)

        end subroutine load_asteroid_database

        subroutine nice_printing()

            print *, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            print *, " "
            print *, "  WELCOME TO THE FORTRAN ASTEROID PREFILTER TOOL (FAP)  "
            print *, " "
            print *, "                  Author: Jack Tyler                     "
            print *, "              Email: jack.tyler@soton.ac.uk              "
            print *, " "
            print *, " This version last updated 2020-02-26                    "
            print *, " "
            print *, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
            print *, " "
            print *, "Parameters specified"
            print *, "~~~~~~~~~~~~~~~~~~~~"
            print *, " "
            print *, "Input data-file: ", target_dataset
            print *, "SPK database", spk_database
            print *, "Target weeks", week_database
            print *,
            print *, "Beginning run "
            print *, "~~~~~~~~~~~~~ "

        end subroutine nice_printing

end module variable_init
