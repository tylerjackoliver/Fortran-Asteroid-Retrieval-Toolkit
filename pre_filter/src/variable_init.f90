module variable_init

    use precision
    use constants
    use problem_parameters

    implicit none
    
    integer, allocatable                        :: database(:)                                           ! SPK ID Database
    integer                                     :: candidate_count = 0                                   ! Number of candidates found
    integer                                     :: itercount = 0                                         ! Main loop sentinel

    real(kind=dp), allocatable                  :: target_state(:,:)                                     ! [r, v] of the target states
    real(kind=dp), allocatable                  :: target_elements(:,:)                                  ! [a, e, i, O, o, M, t] of the target states
    real(kind=dp), allocatable                  :: manifold_data_elements(:,:)                           ! num_targets x 6
    real(kind=dp), allocatable                  :: asteroid_elements(8)                                  ! Orbital elements of the asteroid
    real(kind=dp), allocatable                  :: asteroid_state(6)                                     ! State of the asteroid

    real(kind=dp)                               :: dum                                                   ! Dummy variable

    contains

        subroutine variable_initialisation()

            use precision
            use constants

            implicit none

            integer                             :: iostate
            integer                             :: num_lines   = 0                                             ! Number of lines in the file
            integer                             :: num_targets = 0
            integer                             :: i


            ! Determine number of lines in the target dataset file

            write(*, *, advance='no') "Determining number of lines in the target dataset file..."

            open(unit=99, file=target_dataset)                                                 ! Target OEs
        
            do
        
                read(99, *, iostat=iostate)
        
                if (iostate .ne. 0) then
        
                    exit
        
                else
        
                    num_lines = num_lines + 1
        
                end if

            end do
        
            rewind(99)

            ! Allocate datasets to correct size

            allocate(manifold_data_state(num_lines, 6))
            allocate(manifold_data_elements(num_lines, 8))

            ! Read-in data

            write(*,*, advance='no') "done. Reading dataset into memory..."
        
            do i = 1, num_lines
        
                read(99,*) dum, manifold_data_state(i, :)
                call OSCELT(manifold_data_state(i, :), desireddate, 1.d0, manifold_data_elements(i, :))
        
            end do
        
            close(99)

            write(*,*) "done."
        
            ! Determine the number of asteroids we're targeting

            write(*, *, advance='no') "Determining number of lines in the SPK dataset file..."
        
            open(unit=98, file=spk_database)                                             ! Database file
            
            do
        
                read(98, *, iostat=iostate)
        
                if (iostate .ne. 0) then
        
                    exit
        
                else
        
                    num_targets = num_targets + 1

                end if
        
            end do
        
            rewind(98)

            ! Allocate correct size of database

            allocate(database(num_targets))

            write(*,*, advance='no') "done. Reading dataset into memory..."
        
            do i = 1, num_targets
            
                read(98, *) database(i)
            
            end do
        
            close(98)
        
            open(unit=output_unit, file=output_file)

            ! Load ephemerides

            call FURNSH(planetary_ephemeris)
            call FURNSH(leapseconds_file)

        end subroutine variable_initialisation

        subroutine variable_destruct()

            close(output_unit)
            call UNLOAD(planetary_ephemeris)
            call UNLOAD(leapseconds_file)

            deallocate(database)
            deallocate(manifold_data_elements)
            deallocate(manifold_data_state)

        end subroutine variable_destruct

end module variable_init
