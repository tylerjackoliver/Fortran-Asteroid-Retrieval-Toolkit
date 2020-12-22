module compute_spline

    use bspline_module      ! Uses the spine library
    use problem_parameters  ! Grants access to the dataset

    implicit none

    contains

        ! @brief Interpolates a reference state on the pi/8 planes
        ! @param[in] n_mnfd Value of n_mnfd to interpolate to
        ! @param[in] J Orbit number to interpolate (K in Tyler, J., Wittig, A. 2020: Acta Astronautica, F in Tyler, J., Wittig, A., 2021: AIAA conference)
        subroutine bspline_interpolate(n_mnfd, J, interpolated_state_vector)

            use variable_initialisation

            implicit none

            double precision, intent(in)            :: n_mnfd                           ! Along-orbit direction interpolant
            
            double precision, intent(in)            :: J                                ! Which orbit we're studying (bruting over J)

            double precision, intent(out)           :: interpolated_state_vector(6)     ! Final interpolated structure 

            integer                                 :: bias_n_mnfd                      ! Bias to help with boundary conditions
            integer                                 :: bias_J                           ! Bias to help with boundary conditions
            integer                                 :: nx, ny                           ! Number of abcissae provided
            integer                                 :: kx, ky                           ! Order of spline pieces = 3
            integer                                 :: iknot                            ! Knot sequence flag; set = 0 to automatically choose (here)
            integer                                 :: iflag                            ! Status message; 0 == good
            integer                                 :: J_floor                          ! Floor of t_end
            integer                                 :: J_ceil                           ! Ceiling of t_end
            integer                                 :: n_mnfd_floor                     ! Floor of n_mnfd
            integer                                 :: n_mnfd_ceil                      ! Ceiling of n_mnfd
            integer                                 :: dimension                        ! Dimension of the state vector we're interpolating
            integer                                 :: idx                              ! x-derivative of polynomial to evaluate: set zero to evaluate interpolant
            integer                                 :: idy                              ! As above, but for y
            integer                                 :: inbvx                            ! Initialisation parameter db2val
            integer                                 :: inbvy                            ! Initialisation parameter db2val
            integer                                 :: iloy                             ! Initialisation parameter db2val

            double precision                        :: x_abcissae(4)                    ! Array of x-abcissae
            double precision                        :: y_abcissae(4)                    ! Array of y-abcissae
            double precision                        :: fnc_value(4, 4)                  ! Value of function to use in interpolating
            double precision                        :: tx(7)                            ! Locations of the knots: size (nx + kx)
            double precision                        :: ty(7)                            ! As above, but for the y-coordinate
            double precision                        :: bcoef(4, 4)                      ! Matrix of coefficients of the b-spline interpolant
            double precision                        :: w0(12)                           ! Work array for db2val
            double precision                        :: w1(4)                            ! Work array for db2val

            logical                                 :: extrap = .false.                 ! If extrapolation is OK in db2val

            ! First, initialise problem variables

            nx = 4
            ny = 4
            kx = 3
            ky = 3
            iknot = 0                                                                   ! Decide for me
            iflag = -1                                                                  ! Nothing has been done yet
            idx = 0
            idy = 0

            ! Now, access the correct portions of the dataset for each dimension

            ! Collect the required data
            ! Bias = 0 if we are dealing with 'interior' points
            J_floor = floor(J)
            J_ceil  = ceiling(J)

            n_mnfd_floor = floor(n_mnfd)
            n_mnfd_ceil  = ceiling(n_mnfd)

            ! Check if floor == ceil
            if (J_floor .eq. J_ceil) then

                J_ceil = J_floor + 1

            end if

            if (n_mnfd_floor .eq. n_mnfd_ceil) then

                n_mnfd_ceil = 1 + n_mnfd_floor

            end if

            ! Adjust where we take the grid points on the edges of the grid
            if (J_floor .lt. 2) then ! We need to bias UP

                bias_J = 1

            else if (J_ceil .gt. num_orbits - 1) then ! We need to bias DOWN

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

            do dimension = 2, 7 ! Loop through state vector

                ! Manually override initialisation counters for db2val since new data is computed here
                inbvx = 1
                inbvy = 1
                iloy  = 1

                ! Get function values
                fnc_value(:,:) = dataset( (J_floor+bias_J-1):(J_ceil+bias_J+1), &
                                          (n_mnfd_floor+bias_n_mnfd-1):(n_mnfd_ceil+bias_n_mnfd+1), &
                                            dimension) ! (tend, nmnfd)

                x_abcissae(:) = (/J_floor+bias_J-1, J_floor+bias_J, J_ceil+bias_J, J_ceil+bias_J+1/)
                y_abcissae(:) = (/n_mnfd_floor+bias_n_mnfd-1, n_mnfd_floor+bias_n_mnfd, n_mnfd_ceil+bias_n_mnfd, &
                                n_mnfd_ceil+bias_n_mnfd+1/)

                ! Create the spline
                call db2ink(x_abcissae, nx, y_abcissae, ny, fnc_value, kx, ky, iknot, tx, ty, bcoef, iflag)

                ! If didn't work, exit
                if (iflag .ne. 0) error STOP "Interpolation initialisation failed."

                ! Evaluate the spline in a given dimension
                call db2val(J, n_mnfd, idx, idy, tx, ty, nx, ny, kx, ky, bcoef, interpolated_state_vector(dimension-1), iflag, &
                            inbvx, inbvy, iloy, w1, w0, extrap) 

                ! If didn't work, exit
                if (iflag .ne. 0) error STOP "Interpolation evaluation failed."

            end do ! dimension

        end subroutine bspline_interpolate

        ! @brief Interpolates the time back from the original periodic orbit between grid points
        ! @param[in] J Orbit number to interpolate (K in Tyler, J., Wittig, A. 2020: Acta Astronautica, F in Tyler, J., Wittig, A., 2021: AIAA conference)
        ! @param[in] n_mnfd Value of n_mnfd to interpolate to
        subroutine bspline_interpolate_time(J, n_mnfd, tmani)

            use variable_initialisation

            implicit none
        
            double precision, intent(in)            :: J                                ! Backwards integration time interpolant
            double precision, intent(in)            :: n_mnfd                           ! Along-orbit direction interpolant
            
            double precision, intent(out)           :: tmani                            ! Final interpolated structure 
        
            integer                                 :: bias_n_mnfd                      ! Bias to help with boundary conditions
            integer                                 :: bias_J                           ! Bias to help with boundary conditions
            integer                                 :: nx, ny                           ! Number of abcissae provided
            integer                                 :: kx, ky                           ! Order of spline pieces = 3
            integer                                 :: iknot                            ! Knot sequence flag; set = 0 to automatically choose (here)
            integer                                 :: iflag                            ! Status message; 0 == good
            integer                                 :: J_floor                          ! Floor of J
            integer                                 :: J_ceil                           ! Ceiling of J
            integer                                 :: n_mnfd_floor                     ! Floor of n_mnfd
            integer                                 :: n_mnfd_ceil                      ! Ceiling of n_mnfd
            integer                                 :: dimension                        ! Dimension of the state vector we're interpolating
            integer                                 :: idx                              ! x-derivative of polynomial to evaluate: set zero to evaluate interpolant
            integer                                 :: idy                              ! As above, but for y
            integer                                 :: inbvx                            ! Initialisation parameter db2val
            integer                                 :: inbvy                            ! Initialisation parameter db2val
            integer                                 :: iloy                             ! Initialisation parameter db2val
        
            double precision                        :: x_abcissae(4)                    ! Array of x-abcissae
            double precision                        :: y_abcissae(4)                    ! Array of y-abcissae
            double precision                        :: fnc_value(4, 4)                  ! Value of function to use in interpolating
            double precision                        :: tx(7)                            ! Locations of the knots: size (nx + kx)
            double precision                        :: ty(7)                            ! As above, but for the y-coordinate
            double precision                        :: bcoef(4, 4)                      ! Matrix of coefficients of the b-spline interpolant
            double precision                        :: w0(12)                           ! Work array for db2val
            double precision                        :: w1(4)                            ! Work array for db2val
        
            logical                                 :: extrap = .false.                 ! If extrapolation is OK in db2val
        
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
            J_floor = floor(J)
            J_ceil  = ceiling(J)
        
            n_mnfd_floor = floor(n_mnfd)
            n_mnfd_ceil  = ceiling(n_mnfd)
        
            ! Check if floor == ceil
            if (J_floor .eq. J_ceil) then
        
                J_ceil = J_floor + 1
        
            end if
        
            if (n_mnfd_floor .eq. n_mnfd_ceil) then
        
                n_mnfd_ceil = 1 + n_mnfd_floor
        
            end if
        
            if (J_floor .lt. 2) then ! We need to bias UP
        
                bias_J = 1
        
            else if (J_ceil .gt. num_orbits - 1) then ! We need to bias DOWN
        
                bias_J = -1
        
            else
        
                bias_J = 0
        
            end if
        
            if (n_mnfd_floor .lt. 2) then ! We need to bias UP
        
                bias_n_mnfd = 1
        
            else if (n_mnfd_ceil .gt. n_mnfd_disc - 1) then ! We need to bias DOWN
        
                bias_n_mnfd = -1
        
            else
        
                bias_n_mnfd = 0
        
            end if
        
            do dimension = 1,1 ! Loop through state vector
        
                ! Manually override initialisation counters for db2val since new data is computed here
        
                inbvx = 1
                inbvy = 1
                iloy  = 1
        
                fnc_value(:,:) = dataset( (J_floor+bias_J-1):(J_ceil+bias_J+1), &
                                            (n_mnfd_floor+bias_n_mnfd-1):(n_mnfd_ceil+bias_n_mnfd+1), dimension)
        
                x_abcissae(:) = (/J_floor+bias_J-1, J_floor+bias_J, J_ceil+bias_J, J_ceil+bias_J+1/)
                y_abcissae(:) = (/n_mnfd_floor+bias_n_mnfd-1, n_mnfd_floor+bias_n_mnfd, n_mnfd_ceil+bias_n_mnfd, &
                                n_mnfd_ceil+bias_n_mnfd+1/)
        
                call db2ink(x_abcissae, nx, y_abcissae, ny, fnc_value, kx, ky, iknot, tx, ty, bcoef, iflag)
        
                if (iflag .ne. 0) error STOP "Interpolation initialisation failed."
        
                call db2val(J, n_mnfd, idx, idy, tx, ty, nx, ny, kx, ky, bcoef, tmani, iflag, &
                            inbvx, inbvy, iloy, w1, w0, extrap) 
        
                if (iflag .ne. 0) error STOP "Interpolation evaluation failed."
        
            end do ! dimension
        
        end subroutine bspline_interpolate_time

end module compute_spline
