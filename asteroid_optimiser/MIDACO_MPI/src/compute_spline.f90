module compute_spline

    use precision_kinds     ! Access to working precision types
    use bspline_module      ! Uses the spine library
    use problem_parameters  ! Grants access to the dataset

    implicit none

    contains

        subroutine bspline_interpolate(t_end, n_mnfd, J, interpolated_state_vector)

            implicit none

            real(kind=dp), intent(in)           :: t_end                            ! Backwards integration time interpolant
            real(kind=dp), intent(in)           :: n_mnfd                           ! Along-orbit direction interpolant
            
            integer, intent(in)                 :: J                                ! Which orbit we're studying (bruting over J)

            real(kind=dp), intent(out)          :: interpolated_state_vector(6)     ! Final interpolated structure 

            integer                             :: bias_n_mnfd                      ! Bias to help with boundary conditions
            integer                             :: bias_tend                        ! Bias to help with boundary conditions
            integer                             :: nx, ny                           ! Number of abcissae provided
            integer                             :: kx, ky                           ! Order of spline pieces = 3
            integer                             :: iknot                            ! Knot sequence flag; set = 0 to automatically choose (here)
            integer                             :: iflag                            ! Status message; 0 == good
            integer                             :: t_end_floor                      ! Floor of t_end
            integer                             :: t_end_ceil                       ! Ceiling of t_end
            integer                             :: n_mnfd_floor                     ! Floor of n_mnfd
            integer                             :: n_mnfd_ceil                      ! Ceiling of n_mnfd
            integer                             :: dimension                        ! Dimension of the state vector we're interpolating
            integer                             :: idx                              ! x-derivative of polynomial to evaluate: set zero to evaluate interpolant
            integer                             :: idy                              ! As above, but for y
            integer                             :: inbvx                            ! Initialisation parameter db2val
            integer                             :: inbvy                            ! Initialisation parameter db2val
            integer                             :: iloy                             ! Initialisation parameter db2val

            real(kind=dp)                       :: x_abcissae(4)                    ! Array of x-abcissae
            real(kind=dp)                       :: y_abcissae(4)                    ! Array of y-abcissae
            real(kind=dp)                       :: fnc_value(4, 4)                  ! Value of function to use in interpolating
            real(kind=dp)                       :: tx(7)                            ! Locations of the knots: size (nx + kx)
            real(kind=dp)                       :: ty(7)                            ! As above, but for the y-coordinate
            real(kind=dp)                       :: bcoef(4, 4)                      ! Matrix of coefficients of the b-spline interpolant
            real(kind=dp)                       :: w0(12)                           ! Work array for db2val
            real(kind=dp)                       :: w1(4)                            ! Work array for db2val

            logical                             :: extrap = .false.                 ! If extrapolation is OK in db2val

            ! First, initialise problem variables

            nx = 4
            ny = 4
            kx = 3
            ky = 3
            iknot = 0  ! Decide for me
            iflag = -1 ! Nothing has been done yet
            idx = 0
            idy = 0

            ! Now, access the correct portions of the dataset for each dimension

            ! Collect the required data
            ! Bias = 0 if we are dealing with 'interior' points

            t_end_floor = floor(t_end)
            t_end_ceil  = ceiling(t_end)

            n_mnfd_floor = floor(n_mnfd)
            n_mnfd_ceil  = ceiling(n_mnfd)

            ! Check if floor == ceil

            if (t_end_floor .eq. t_end_ceil) then

                t_end_ceil = t_end_floor + 1

            end if

            if (n_mnfd_floor .eq. n_mnfd_ceil) then

                n_mnfd_ceil = 1 + n_mnfd_floor

            end if

            if (t_end_floor .lt. 2) then ! We need to bias UP

                bias_tend = 1

            else if (t_end_ceil .gt. 99) then ! We need to bias DOWN

                bias_tend = -1

            else

                bias_tend = 0

            end if

            if (n_mnfd_floor .lt. 2) then ! We need to bias UP

                bias_n_mnfd = 1

            else if (n_mnfd_ceil .gt. 99) then ! We need to bias DOWN

                bias_n_mnfd = -1

            else

                bias_n_mnfd = 0

            end if

            do dimension = 1,6 ! Loop through state vector

                ! Manually override initialisation counters for db2val since new data is computed here

                inbvx = 1
                inbvy = 1
                iloy  = 1

                fnc_value(:,:) = dataset( (t_end_floor+bias_tend-1):(t_end_ceil+bias_tend+1), &
                                          (n_mnfd_floor+bias_n_mnfd-1):(n_mnfd_ceil+bias_n_mnfd+1), &
                                            J, dimension) ! (tend, nmnfd)

                x_abcissae(:) = (/t_end_floor+bias_tend-1, t_end_floor+bias_tend, t_end_ceil+bias_tend, t_end_ceil+bias_tend+1/)
                y_abcissae(:) = (/n_mnfd_floor+bias_n_mnfd-1, n_mnfd_floor+bias_n_mnfd, n_mnfd_ceil+bias_n_mnfd, &
                                n_mnfd_ceil+bias_n_mnfd+1/)

                call db2ink(x_abcissae, nx, y_abcissae, ny, fnc_value, kx, ky, iknot, tx, ty, bcoef, iflag)

                if (iflag .ne. 0) error STOP "Interpolation initialisation failed."

                call db2val(t_end, n_mnfd, idx, idy, tx, ty, nx, ny, kx, ky, bcoef, interpolated_state_vector(dimension), iflag, &
                            inbvx, inbvy, iloy, w1, w0, extrap) 

                if (iflag .ne. 0) error STOP "Interpolation evaluation failed."

            end do ! dimension

        end subroutine bspline_interpolate

        subroutine b_spline_interpolate_perturbed_conds(n_mnfd, J, interpolated_state_vector)

            implicit none

            real(kind=dp), intent(in)           :: n_mnfd                           ! Along-orbit direction interpolant
            
            integer, intent(in)                 :: J                                ! Which orbit we're studying (bruting over J)

            real(kind=dp)                       :: temp_state_vector(7)             ! Holds time as well
            real(kind=dp), intent(out)          :: interpolated_state_vector(6)     ! Final interpolated structure 

            integer                             :: bias_n_mnfd                      ! Bias to help with boundary conditions
            integer                             :: nx                               ! Number of abcissae provided
            integer                             :: kx                               ! Order of spline pieces = 3
            integer                             :: iknot                            ! Knot sequence flag; set = 0 to automatically choose (here)
            integer                             :: iflag                            ! Status message; 0 == good
            integer                             :: n_mnfd_floor                     ! Floor of n_mnfd
            integer                             :: n_mnfd_ceil                      ! Ceiling of n_mnfd
            integer                             :: dimension                        ! Dimension of the state vector we're interpolating
            integer                             :: idx                              ! x-derivative of polynomial to evaluate: set zero to evaluate interpolant
            integer                             :: inbvx                            ! Initialisation parameter db2val
            integer                             :: iloy                             ! Initialisation parameter db2val

            real(kind=dp)                       :: x_abcissae(4)                    ! Array of x-abcissae
            real(kind=dp)                       :: fnc_value(4)                     ! Value of function to use in interpolating
            real(kind=dp)                       :: tx(7)                            ! Locations of the knots: size (nx + kx)
            real(kind=dp)                       :: bcoef(4)                         ! Matrix of coefficients of the b-spline interpolant
            real(kind=dp)                       :: w0(12)                           ! Work array for db2val
            real(kind=dp)                       :: w1(4)                            ! Work array for db2val

            logical                             :: extrap = .false.                 ! If extrapolation is OK in db2val

            ! First, initialise problem variables

            nx = 4
            kx = 3
            iknot = 0  ! Decide for me
            iflag = -1 ! Nothing has been done yet
            idx = 0

            ! Now, access the correct portions of the dataset for each dimension

            ! Collect the required data
            ! Bias = 0 if we are dealing with 'interior' points

            n_mnfd_floor = floor(n_mnfd)
            n_mnfd_ceil  = ceiling(n_mnfd)

            ! Check if floor == ceil

            if (n_mnfd_floor .eq. n_mnfd_ceil) then

                n_mnfd_ceil = 1 + n_mnfd_floor

            end if

            if (n_mnfd_floor .lt. 2) then ! We need to bias UP

                bias_n_mnfd = 1

            else if (n_mnfd_ceil .gt. 99) then ! We need to bias DOWN

                bias_n_mnfd = -1

            else

                bias_n_mnfd = 0

            end if

            do dimension = 1,7 ! Loop through state vector

                ! Manually override initialisation counters for db2val since new data is computed here

                inbvx = 1
                inbvy = 1
                iloy  = 1

                fnc_value(:,:) = perturbed_conds_dataset( (n_mnfd_floor+bias_n_mnfd-1):(n_mnfd_ceil+bias_n_mnfd+1), &
                                            J, dimension) ! (tend, nmnfd)

                x_abcissae(:) = (/n_mnfd_floor+bias_n_mnfd-1, n_mnfd_floor+bias_n_mnfd, n_mnfd_ceil+bias_n_mnfd, &
                                n_mnfd_ceil+bias_n_mnfd+1/)

                call db1ink(x_abcissae, nx, y_abcissae, ny, fnc_value, kx, ky, iknot, tx, ty, bcoef, iflag)

                if (iflag .ne. 0) error STOP "Interpolation initialisation failed."

                call db1val(n_mnfd, idx, tx, nx, kx, bcoef, temp_state_vector(dimension), iflag, &
                            inbvx, w0, extrap) 

                if (iflag .ne. 0) error STOP "Interpolation evaluation failed."

            end do ! dimension

            tmani = temp_state_vector(1)
            interpolated_state_vector = temp_state_vector(2:7)

        end subroutine b_spline_interpolate_perturbed_conds

end module compute_spline
