!------------------------------------------------------------------------------
! Fortran Asteroid Retrieval Tool (FART) v1.0: state_determination module
!------------------------------------------------------------------------------
!
! MODULE: State determination
!
!> @author
!> Jack Tyler, University of Southampton
!
! DESCRIPTION: 
!> Contains routines to construct interpolants for the optimisation routines.
!
! REVISION HISTORY:
! 01 Mar 2018 - Initial Version
! 01 Jul 2020 - Refactoring; add Doxygen support
!------------------------------------------------------------------------------


module compute_spline

    use bspline_module              ! Use the BSPLINE-FORTRAN library
    use problem_parameters          ! Grants access to the dataset and other problem constants
    use variable_initialisation     ! Grants access to problem constants

    implicit none

    public :: bspline_interpolate, b_spline_interpolate_perturbed_conds

    contains

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Interpolates between two non-integer values of n_mnfd and J on the +-pi/8
        !! plane.
        !
        ! REVISION HISTORY:
        ! 01 Feb 2020 - Initial version
        ! 04 July 2020 - Refactored; Doxygen support
        !
        !> @param[in] n_mnfd, J
        !> @param[out] interpolated_state_vector
        !--------------------------------------------------------------------------- 

        subroutine bspline_interpolate(n_mnfd, J, interpolated_state_vector)

            implicit none

            double precision, intent(in)        :: n_mnfd                           ! Intra-orbit interpolant
            double precision, intent(in)        :: J                                ! Inter-orbit interpolant

            double precision, intent(out)       :: interpolated_state_vector(6)     ! Final interpolated state vector

            integer                             :: bias_n_mnfd                      ! Bias to assist boundary conditions
            integer                             :: bias_J                           ! Bias to assist boundary conditions
            integer                             :: nx, ny                           ! Number of abcissae provided
            integer                             :: kx, ky                           ! Order of spline pieces = 3
            integer                             :: iknot                            ! Knot sequence flag; set = 0 to automatically choose (here)
            integer                             :: iflag                            ! Status message; 0 == good
            integer                             :: J_floor                          ! Floor of J (handling BCs)
            integer                             :: J_ceil                           ! Ceiling of J (handling BCs)
            integer                             :: n_mnfd_floor                     ! Floor of n_mnfd (handling BCs)
            integer                             :: n_mnfd_ceil                      ! Ceiling of n_mnfd (handling BCs)
            integer                             :: dimension                        ! Dimension of the state vector we're interpolating
            integer                             :: idx                              ! x-derivative of polynomial to evaluate: set zero to evaluate interpolant
            integer                             :: idy                              ! As above, but for y
            integer                             :: inbvx                            ! Initialisation parameter db2val
            integer                             :: inbvy                            ! Initialisation parameter db2val
            integer                             :: iloy                             ! Initialisation parameter db2val

            double precision                    :: x_abcissae(4)                    ! Array of x-abcissae
            double precision                    :: y_abcissae(4)                    ! Array of y-abcissae
            double precision                    :: fnc_value(4, 4)                  ! Value of function to use in interpolating
            double precision                    :: tx(7)                            ! Locations of the knots: size (nx + kx)
            double precision                    :: ty(7)                            ! As above, but for the y-coordinate
            double precision                    :: bcoef(4, 4)                      ! Matrix of coefficients of the b-spline interpolant
            double precision                    :: w0(12)                           ! Work array for db2val
            double precision                    :: w1(4)                            ! Work array for db2val

            logical                             :: extrap = .false.                 ! If extrapolation is OK in db2val

            ! First, initialise problem variables

            nx = 4
            ny = 4
            kx = 3
            ky = 3
            iknot = 0  ! Decide location of the knots for me
            iflag = -1 ! Nothing has been done yet
            idx = 0
            idy = 0

            !
            ! Use the floors/ceils of J, n_mnfd to determine upper and lower ranges of the dataset needed
            ! to correctly determine the interpolant.
            !
            ! By default, 4 points are used. Thus, we have: J range=(floor(J) - 1, floor(J), ceil(J), ceil(J)+1), and
            ! similarly for n_mnfd
            !

            J_floor = floor(J)
            J_ceil  = ceiling(J)

            n_mnfd_floor = floor(n_mnfd)
            n_mnfd_ceil  = ceiling(n_mnfd)

            !
            ! Check if floor == ceil; if so, add one to ceil to adjust
            !

            if (J_floor .eq. J_ceil) then

                J_ceil = J_floor + 1

            end if

            if (n_mnfd_floor .eq. n_mnfd_ceil) then

                n_mnfd_ceil = 1 + n_mnfd_floor

            end if

            !
            ! Now, set biases. These adjust the range of data 'up' and 'down' in each variable to handle BCs.
            ! For example, if we want orbit J = 1.1, then we cannot obtain floor(J)-1. In this case, we use J=(1,2,3,4),
            ! so we add a bias of 1.
            !

            if (J_floor .lt. 2) then ! We need to bias UP

                bias_J = 1

            else if (J_ceil .gt. num_orbits-1) then ! We need to bias DOWN

                bias_J = -1

            else

                bias_J = 0

            end if

            if (n_mnfd_floor .lt. 2) then ! We need to bias UP

                bias_n_mnfd = 1

            else if (n_mnfd_ceil .gt. n_mnfd_disc-1) then ! We need to bias DOWN

                bias_n_mnfd = -1

            else

                bias_n_mnfd = 0

            end if

            !
            ! Now construct interpolants for each dimension in turn, and add this to the interpolated state vector
            !

            do dimension = 2, 7 ! Loop through state vector; start at 2 as position 1 is time

                ! Manually override initialisation counters for db2val since new data is computed here

                inbvx = 1
                inbvy = 1
                iloy  = 1

                !
                ! Now populate abcissae arrays
                !

                x_abcissae(:) = (/J_floor+bias_J-1, J_floor+bias_J, J_ceil+bias_J, J_ceil+bias_J+1/)
                y_abcissae(:) = (/n_mnfd_floor+bias_n_mnfd-1, n_mnfd_floor+bias_n_mnfd, n_mnfd_ceil+bias_n_mnfd, &
                                n_mnfd_ceil+bias_n_mnfd+1/)

                !
                ! Obtain data points using slicing convention
                !

                fnc_value(:,:) = dataset( x_abcissae(1):x_abcissae(4), &
                                          y_abcissae(1):y_abcissae(4), &
                                          dimension )

                !
                ! Initialise interpolant construction
                !

                call db2ink(x_abcissae, nx, y_abcissae, ny, fnc_value, kx, ky, iknot, tx, ty, bcoef, iflag)

                if (iflag .ne. 0) error STOP "Interpolation initialisation failed." ! If error, stop

                !
                ! Compute interpolant
                !

                call db2val(J, n_mnfd, idx, idy, tx, ty, nx, ny, kx, ky, bcoef, interpolated_state_vector(dimension-1), iflag, &
                            inbvx, inbvy, iloy, w1, w0, extrap) 

                if (iflag .ne. 0) error STOP "Interpolation evaluation failed." ! If error, stop

            end do ! dimension

        end subroutine bspline_interpolate

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of routine. 
        !> @brief
        !> Interpolates between two non-integer values of n_mnfd and J in the initial
        !! conditions for orbits.
        !
        ! REVISION HISTORY:
        ! 01 Feb 2020 - Initial version
        ! 04 July 2020 - Refactored; Doxygen support
        !
        !> @param[in] n_mnfd, J
        !> @param[out] interpolated_state_vector
        !--------------------------------------------------------------------------- 

        subroutine b_spline_interpolate_perturbed_conds(n_mnfd, J, interpolated_state_vector)
            
            use variable_initialisation

            implicit none

            double precision, intent(in)        :: n_mnfd                           ! Intra-orbit interpolant
            double precision, intent(in)        :: J                                ! Inter-orbit interpolant

            double precision, intent(out)       :: interpolated_state_vector(6)     ! Final interpolated state vector

            integer                             :: bias_n_mnfd                      ! Bias to assist boundary conditions
            integer                             :: bias_J                           ! Bias to assist boundary conditions
            integer                             :: nx, ny                           ! Number of abcissae provided
            integer                             :: kx, ky                           ! Order of spline pieces = 3
            integer                             :: iknot                            ! Knot sequence flag; set = 0 to automatically choose (here)
            integer                             :: iflag                            ! Status message; 0 == good
            integer                             :: J_floor                          ! Floor of J (handling BCs)
            integer                             :: J_ceil                           ! Ceiling of J (handling BCs)
            integer                             :: n_mnfd_floor                     ! Floor of n_mnfd (handling BCs)
            integer                             :: n_mnfd_ceil                      ! Ceiling of n_mnfd (handling BCs)
            integer                             :: dimension                        ! Dimension of the state vector we're interpolating
            integer                             :: idx                              ! x-derivative of polynomial to evaluate: set zero to evaluate interpolant
            integer                             :: idy                              ! As above, but for y
            integer                             :: inbvx                            ! Initialisation parameter db2val
            integer                             :: inbvy                            ! Initialisation parameter db2val
            integer                             :: iloy                             ! Initialisation parameter db2val

            double precision                    :: x_abcissae(4)                    ! Array of x-abcissae
            double precision                    :: y_abcissae(4)                    ! Array of y-abcissae
            double precision                    :: fnc_value(4, 4)                  ! Value of function to use in interpolating
            double precision                    :: tx(7)                            ! Locations of the knots: size (nx + kx)
            double precision                    :: ty(7)                            ! As above, but for the y-coordinate
            double precision                    :: bcoef(4, 4)                      ! Matrix of coefficients of the b-spline interpolant
            double precision                    :: w0(12)                           ! Work array for db2val
            double precision                    :: w1(4)                            ! Work array for db2val

            logical                             :: extrap = .false.                 ! If extrapolation is OK in db2val

            ! First, initialise problem variables

            nx = 4
            ny = 4
            kx = 3
            ky = 3
            iknot = 0  ! Decide location of the knots for me
            iflag = -1 ! Nothing has been done yet
            idx = 0
            idy = 0

            !
            ! Use the floors/ceils of J, n_mnfd to determine upper and lower ranges of the dataset needed
            ! to correctly determine the interpolant.
            !
            ! By default, 4 points are used. Thus, we have: J range=(floor(J) - 1, floor(J), ceil(J), ceil(J)+1), and
            ! similarly for n_mnfd
            !

            J_floor = floor(J)
            J_ceil  = ceiling(J)

            n_mnfd_floor = floor(n_mnfd)
            n_mnfd_ceil  = ceiling(n_mnfd)

            !
            ! Check if floor == ceil; if so, add one to ceil to adjust
            !

            if (J_floor .eq. J_ceil) then

                J_ceil = J_floor + 1

            end if

            if (n_mnfd_floor .eq. n_mnfd_ceil) then

                n_mnfd_ceil = 1 + n_mnfd_floor

            end if

            !
            ! Now, set biases. These adjust the range of data 'up' and 'down' in each variable to handle BCs.
            ! For example, if we want orbit J = 1.1, then we cannot obtain floor(J)-1. In this case, we use J=(1,2,3,4),
            ! so we add a bias of 1.
            !

            if (J_floor .lt. 2) then ! We need to bias UP

                bias_J = 1

            else if (J_ceil .gt. num_orbits-1) then ! We need to bias DOWN

                bias_J = -1

            else

                bias_J = 0

            end if

            if (n_mnfd_floor .lt. 2) then ! We need to bias UP

                bias_n_mnfd = 1

            else if (n_mnfd_ceil .gt. n_mnfd_disc-1) then ! We need to bias DOWN

                bias_n_mnfd = -1

            else

                bias_n_mnfd = 0

            end if

            do dimension = 1, 6 ! Loop through state vector

                ! Manually override initialisation counters for db2val since new data is computed here

                inbvx = 1
                inbvy = 1
                iloy  = 1

                x_abcissae(:) = (/n_mnfd_floor+bias_n_mnfd-1, n_mnfd_floor+bias_n_mnfd, n_mnfd_ceil+bias_n_mnfd, &
                n_mnfd_ceil+bias_n_mnfd+1/)

                y_abcissae(:) = (/J_floor+bias_J-1, J_floor+bias_J, J_ceil+bias_J, J_ceil+bias_J+1/)

                fnc_value(:,:) = perturbed_conds_dataset(   x_abcissae(1):x_abcissae(4), &
                                                            y_abcissae(1):y_abcissae(4), &
                                                            dimension )

                !
                ! Initialise interpolant
                !

                call db2ink(x_abcissae, nx, y_abcissae, ny, fnc_value, kx, ky, iknot, tx, ty, bcoef, iflag)

                if (iflag .ne. 0) error STOP "Interpolation initialisation failed." ! If error, stop

                !
                ! Evaluate interpolant
                !

                call db2val(n_mnfd, J, idx, idy, tx, ty, nx, ny, kx, ky, bcoef, interpolated_state_vector(dimension), iflag, &
                            inbvx, inbvy, iloy, w1, w0, extrap) 

                if (iflag .ne. 0) error STOP "Interpolation evaluation failed." ! If error

            end do ! dimension

        end subroutine b_spline_interpolate_perturbed_conds

end module compute_spline
