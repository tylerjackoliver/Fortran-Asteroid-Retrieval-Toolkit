module local_optimisation

     
    use constants
    use problem_parameters
    use variable_initialisation
    use midaco_interface
    use ancillary_data
    use sorting_routines

    implicit none

    ! Global arrays - used for capturing the best candidates for global optimisation. Double precision array of n x 6 (dV + optimisation variables)

    double precision, allocatable :: best_optimisation_candidates(:,:) ! n x 6
    double precision, allocatable :: lower_bounds(:,:)
    double precision, allocatable :: upper_bounds(:,:)

    contains

        subroutine initialise_local_optim(number_of_local_solutions)

            ! Inputs

            integer, intent(in)     :: number_of_local_solutions ! How many solutions in the Pareto front do we want to investigate in the local solver?

            integer                 :: solution_number

            !
            ! Subroutine should wrap all required functions to end global processing, and then begin
            ! computing a local approximation
            !
            ! Therefore:    - Should get the pareto front for the current global optimisation
            !               - Should then extract the (10?) most promising points from the Pareto front to move forward into local optimisation
            !               - Reset MIDACO and provide new best-case runs
            !               - *Then* do the intermediate variable destruction 
            !

            ! Get the pareto front for the current global optimisation

            call get_pareto_front()

            ! Get the best candidates for local optimisation

            call get_best_points(number_of_local_solutions)

            ! Initialise the upper and lower bounds

            call initialise_bounds(number_of_local_solutions)

        end subroutine initialise_local_optim

        subroutine initialise_bounds(number_of_local_solutions)

            integer, intent(in) :: number_of_local_solutions

            integer             :: solution_number
            integer             :: dim
            double precision    :: tmp_swap

            double precision    :: reference_lower_bound(5)
            double precision    :: reference_upper_bound(5)

            call STR2ET('2025 Jan 01 00:00', reference_lower_bound(1))
            reference_lower_bound(2) = 1.d0 * 86400.
            reference_lower_bound(3) = -250.d0
            reference_lower_bound(4) = 2
            reference_lower_bound(5) = 2

            call STR2ET('2099 Dec 31 00:00', reference_upper_bound(1))
            reference_upper_bound(2) = 1400.d0 * 86400
            reference_upper_bound(3) = 0.d0                                               ! t_end - should be the *most negative* time of all pi/8 planes
            reference_upper_bound(4) = n_mnfd_disc - 1
            reference_upper_bound(5) = num_orbits - 1

            allocate(lower_bounds(number_of_local_solutions, 5))
            allocate(upper_bounds(number_of_local_solutions, 5))

            do solution_number = 1, number_of_local_solutions

                lower_bounds(solution_number, :) = 0.95 * best_optimisation_candidates(solution_number, 2:6)
                upper_bounds(solution_number, :) = 1.05 * best_optimisation_candidates(solution_number, 2:6)

                tmp_swap = lower_bounds(solution_number, 3)

                lower_bounds(solution_number, 3) = upper_bounds(solution_number, 3) ! As 3 is negative
                upper_bounds(solution_number, 3) = tmp_swap

                ! Fix occasions where 1.2 or 0.8 x the solution would infringe on lower or upper bounds from the dataset

                do dim = 1, 5

                    if (lower_bounds(solution_number, dim) .lt. reference_lower_bound(dim)) then
                        
                        lower_bounds(solution_number, dim) = reference_lower_bound(dim)
                    
                    end if
                    
                    if (upper_bounds(solution_number, dim) .gt. reference_upper_bound(dim)) then
                        
                        upper_bounds(solution_number, dim) = reference_upper_bound(dim)

                    end if

                end do

            end do

        end subroutine initialise_bounds

        subroutine initialise_local_variables(run_number)

            implicit none

            integer :: temp
            integer, intent(in) :: run_number

            character(1023) :: run_number_str

            write(run_number_str, '(I0)') run_number
            ! call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
            ! error stop

            ! Reset MIDACO

            temp = reset_midaco_run()

            ! Close the current pareto front file - open a new one

            close(iunit)

            open(unit=iunit, file='../refined_solutions/paretoFront_'//trim(targ_can)//run_number_str)!//trim(run_number_str))

            ! Optimiser bounds - set crudely as +- 20% from current solution

            XL = lower_bounds(run_number, :)
            XU = upper_bounds(run_number, :)

            ! Starting point, XOPT

            XOPT = (XL + XU)/2.d0

            ! Maximum function evaluations

            max_eval = 99999999
            optim_flag = 0
            optim_stop = 0

            local_time = max_time_local

            ! Printing options

            print_eval = 1000
            save_to_file = 0 ! Save solution to text files

            ! Choose MIDACO parameters (FOR ADVANCED USERS)

            PARAM( 1) = 0.0D0           ! ACCURACY
            PARAM( 2) = 0.0D0           ! SEED
            PARAM( 3) = 0.0D0           ! FSTOP
            PARAM( 4) = 0.0D0           ! ALGOSTOP
            PARAM( 5) = 0.0D0           ! EVALSTOP
            PARAM( 6) = 100000.D0       ! FOCUS
            PARAM( 7) = 1000.0D0        ! ANTS
            PARAM( 8) = 50.0D0          ! KERNEL
            PARAM( 9) = 0.0D0           ! ORACLE
            PARAM(10) = 1000D0          ! PARETOMAX
            PARAM(11) = 0.001D0         ! EPSILON  
            PARAM(12) = -1.000D0        ! BALANCE => focus only on first objective function (DeltaV, not tt)
            PARAM(13) = 0.0D0           ! CHARACTER 

            iw = 0
            rw = 0.d0
            pf = 0.d0

            temp = reset_midaco_run()

        end subroutine initialise_local_variables

        subroutine destroy_local_variables()

            close(iunit)

        end subroutine destroy_local_variables

        subroutine get_best_points(number_of_local_solutions)

            ! Inputs

            integer, intent(in)     :: number_of_local_solutions ! How many solutions in the Pareto front do we want to investigate in the local solver?

            ! Local arrays

            double precision, allocatable :: pareto_points(:,:)   ! Temporary 'work' array - holds all the Pareto points found ready for sorting
            
            integer, allocatable    :: permutation_index(:) ! Array that holds the permutation index relating to the best $n$ Pareto points

            integer                 :: solution_number
            integer                 :: variable_iter
            integer                 :: number_of_pareto_solutions
            integer                 :: idx
            integer                 :: PARETOMAX

            ! Allocate storage space for the optimisation candidates based on the total number of solutions desired

            allocate(best_optimisation_candidates(number_of_local_solutions, 6))

            ! Get the local solutions

            number_of_pareto_solutions = PF(1)  ! First element stores the number of Pareto solutions
            PARETOMAX = PARAM(10)               ! Defined from PARETO configuration arrays

            ! Allocate storage space for all of the Pareto candidates based on the total number of solutions that exist

            allocate(pareto_points(number_of_pareto_solutions, 6))
            allocate(permutation_index(number_of_local_solutions))

            ! Now we know the number of Pareto solution, we can get the solutions

            do solution_number = 1, number_of_pareto_solutions

                pareto_points(solution_number, 1) = PF(2 + O*(solution_number-1)) ! O is defined globally - number of objectives
                pareto_points(solution_number, 2:6) = PF( (2 + O * PARETOMAX + M * PARETOMAX + &
                                                            n*(solution_number-1)) : (2 + O * PARETOMAX + M * PARETOMAX + &
                                                            n*(solution_number-1) + 4) ) ! See MIDACO solver manual for definitions of each term + structure of PF

            end do

            ! Sort the solutions array - get permutations

            call sort_array(number_of_local_solutions, pareto_points(:, 1), permutation_index)

            ! Now it has been sorted, isolate the top 5!

            do solution_number = 1, number_of_local_solutions

                idx = permutation_index(solution_number)

                ! Get this entry in pareto_points; add to the best_candidates array

                best_optimisation_candidates(solution_number, :) = pareto_points(idx, :)


            end do

        end subroutine get_best_points

end module local_optimisation
    