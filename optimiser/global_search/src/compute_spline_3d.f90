module compute_spline

          ! Access to working precision types
    use bspline_module      ! Uses the spine library
    use problem_parameters  ! Grants access to the dataset

    implicit none

    contains

        subroutine bspline_interpolate(t_end, n_mnfd, J, interpolated_state_vector)

            use variable_initialisation

            implicit none

            double precision, intent(in)           :: t_end                            ! Backwards integration time interpolant
            double precision, intent(in)           :: n_mnfd                           ! Along-orbit direction interpolant
            
            double precision, intent(in)           :: J                                ! Which orbit we're studying (bruting over J)

            double precision, intent(out)          :: interpolated_state_vector(6)     ! Final interpolated structure 

            integer                             :: bias_n_mnfd                      ! Bias to help with boundary conditions
            integer                             :: bias_tend                        ! Bias to help with boundary conditions
            integer                             :: bias_J                           ! Bias to help with boundary conditions
            integer                             :: nx, ny, nz                       ! Number of abcissae provided
            integer                             :: kx, ky, kz                       ! Order of spline pieces = 3
            integer                             :: iknot                            ! Knot sequence flag; set = 0 to automatically choose (here)
            integer                             :: iflag                            ! Status message; 0 == good
            integer                             :: t_end_floor                      ! Floor of t_end
            integer                             :: t_end_ceil                       ! Ceiling of t_end
            integer                             :: n_mnfd_floor                     ! Floor of n_mnfd
            integer                             :: n_mnfd_ceil                      ! Ceiling of n_mnfd
            integer                             :: J_floor                          ! Floor of the desired orbit
            integer                             :: J_ceil                           ! Ceiling of the desired orbit
            integer                             :: dimension                        ! Dimension of the state vector we're interpolating
            integer                             :: idx                              ! x-derivative of polynomial to evaluate: set zero to evaluate interpolant
            integer                             :: idy                              ! As above, but for y
            integer                             :: idz                              ! As above, but for z
            integer                             :: inbvx                            ! Initialisation parameter db3val
            integer                             :: inbvy                            ! Initialisation parameter db3val
            integer                             :: inbvz                            ! Initialisation parameter db3val
            integer                             :: iloy                             ! Initialisation parameter db3val
            integer                             :: iloz                             ! Initialisation parameter db3val

            double precision                       :: x_abcissae(4)                    ! Array of x-abcissae
            double precision                       :: y_abcissae(4)                    ! Array of y-abcissae
            double precision                       :: z_abcissae(4)                    ! Array of z-abcissae
            double precision                       :: fnc_value(4, 4, 4)               ! Value of function to use in interpolating
            double precision                       :: tx(7)                            ! Locations of the knots: size (nx + kx)
            double precision                       :: ty(7)                            ! As above, but for the y-coordinate
            double precision                       :: tz(7)                            ! As above, but for the z-coordinate
            double precision                       :: bcoef(4, 4, 4)                   ! Matrix of coefficients of the b-spline interpolant
            double precision                       :: w0(12)                           ! Work array for db2val
            double precision                       :: w1(4)                            ! Work array for db2val
            double precision                       :: w2(3, 3)                         ! Work array db3val

            logical                             :: extrap = .false.                 ! If extrapolation is OK in db2val

            ! First, initialise problem variables

            nx = 4
            ny = 4
            nz = 4
            kx = 3
            ky = 3
            kz = 3
            iknot = 0  ! Decide for me
            iflag = -1 ! Nothing has been done yet
            idx = 0
            idy = 0
            idz = 0

            ! Now, access the correct portions of the dataset for each dimension

            ! Collect the required data
            ! Bias = 0 if we are dealing with 'interior' points

            t_end_floor = floor(t_end)
            t_end_ceil  = ceiling(t_end)

            n_mnfd_floor = floor(n_mnfd)
            n_mnfd_ceil  = ceiling(n_mnfd)

            J_floor     = floor(J) 
            J_ceil      = ceiling(J)

            ! Check if floor == ceil

            if (t_end_floor .eq. t_end_ceil) then

                t_end_ceil = t_end_floor + 1

            end if

            if (n_mnfd_floor .eq. n_mnfd_ceil) then

                n_mnfd_ceil = 1 + n_mnfd_floor

            end if

            if (J_floor .eq. J_ceil) then

                J_ceil = J_floor + 1

            end if

            if (t_end_floor .lt. 2) then ! We need to bias UP

                bias_tend = 1

            else if (t_end_ceil .gt. t_end_disc-1) then ! We need to bias DOWN

                bias_tend = -1

            else

                bias_tend = 0

            end if

            if (n_mnfd_floor .lt. 2) then ! We need to bias UP

                bias_n_mnfd = 1

            else if (n_mnfd_ceil .gt. n_mnfd_disc - 1) then ! We need to bias DOWN

                bias_n_mnfd = -1

            else

                bias_n_mnfd = 0

            end if

            if (J_floor .lt. 2) then

                bias_J = 1

            else if (J_ceil .gt. num_orbits - 1) then

                bias_J = -1

            else

                bias_J = 0

            end if

            do dimension = 2, 7 ! Loop through state vector

                ! Manually override initialisation counters for db2val since new data is computed here

                inbvx = 1
                inbvy = 1
                inbvz = 1
                iloy  = 1
                iloz  = 1

                fnc_value(:,:, :) = dataset( (t_end_floor+bias_tend-1):(t_end_ceil+bias_tend+1), &
                                          (n_mnfd_floor+bias_n_mnfd-1):(n_mnfd_ceil+bias_n_mnfd+1), &
                                            (J_floor+bias_J-1):(J_ceil+bias_J+1), dimension) ! (tend, nmnfd)

                x_abcissae(:) = (/t_end_floor+bias_tend-1, t_end_floor+bias_tend, t_end_ceil+bias_tend, t_end_ceil+bias_tend+1/)
                y_abcissae(:) = (/n_mnfd_floor+bias_n_mnfd-1, n_mnfd_floor+bias_n_mnfd, n_mnfd_ceil+bias_n_mnfd, &
                                n_mnfd_ceil+bias_n_mnfd+1/)
                z_abcissae(:) = (/J_floor+bias_J-1, J_floor+bias_J, J_ceil+bias_J, J_ceil+bias_J+1/)

                call db3ink(x_abcissae, nx, y_abcissae, ny, z_abcissae, nz, fnc_value, kx, ky, kz, iknot, tx, ty, tz, bcoef, &
                            iflag)

                if (iflag .ne. 0) then
                    
                    print *, "iflag at failure", iflag
                    error STOP "Interpolation initialisation failed."

                end if

                call db3val(t_end, n_mnfd, J, idx, idy, idz, tx, ty, tz, nx, ny, nz, kx, ky, kz, bcoef, &
                            interpolated_state_vector(dimension-1), iflag, inbvx, inbvy, inbvz, iloy, iloz, w2, w1, w0, extrap) 

                if (iflag .ne. 0) then
                    
                    print *, "iflag at failure", iflag
                    error STOP "Interpolation evaluation failed."

                end if

            end do ! dimension

        end subroutine bspline_interpolate


        subroutine bspline_interpolate_step(t_end, n_mnfd, J, interpolated_state_vector)

            use variable_initialisation

            implicit none

            double precision, intent(in)           :: t_end                            ! Backwards integration time interpolant
            double precision, intent(in)           :: n_mnfd                           ! Along-orbit direction interpolant
            
            double precision, intent(in)           :: J                                ! Which orbit we're studying (bruting over J)

            double precision, intent(out)          :: interpolated_state_vector(6)     ! Final interpolated structure 

            integer                             :: bias_n_mnfd                      ! Bias to help with boundary conditions
            integer                             :: bias_tend                        ! Bias to help with boundary conditions
            integer                             :: bias_J                           ! Bias to help with boundary conditions
            integer                             :: nx, ny, nz                       ! Number of abcissae provided
            integer                             :: kx, ky, kz                       ! Order of spline pieces = 3
            integer                             :: iknot                            ! Knot sequence flag; set = 0 to automatically choose (here)
            integer                             :: iflag                            ! Status message; 0 == good
            integer                             :: t_end_floor                      ! Floor of t_end
            integer                             :: t_end_ceil                       ! Ceiling of t_end
            integer                             :: n_mnfd_floor                     ! Floor of n_mnfd
            integer                             :: n_mnfd_ceil                      ! Ceiling of n_mnfd
            integer                             :: J_floor                          ! Floor of the desired orbit
            integer                             :: J_ceil                           ! Ceiling of the desired orbit
            integer                             :: dimension                        ! Dimension of the state vector we're interpolating
            integer                             :: idx                              ! x-derivative of polynomial to evaluate: set zero to evaluate interpolant
            integer                             :: idy                              ! As above, but for y
            integer                             :: idz                              ! As above, but for z
            integer                             :: inbvx                            ! Initialisation parameter db3val
            integer                             :: inbvy                            ! Initialisation parameter db3val
            integer                             :: inbvz                            ! Initialisation parameter db3val
            integer                             :: iloy                             ! Initialisation parameter db3val
            integer                             :: iloz                             ! Initialisation parameter db3val
            integer                             :: ii, jj, kk

            double precision                       :: x_abcissae(4)                    ! Array of x-abcissae
            double precision                       :: x_abcissae_step(4)               ! Array of x-abcissae for the time-steps
            double precision                       :: y_abcissae(4)                    ! Array of y-abcissae
            double precision                       :: z_abcissae(4)                    ! Array of z-abcissae
            double precision                       :: fnc_value(4, 4, 4)               ! Value of function to use in interpolating
            double precision                       :: tx(7)                            ! Locations of the knots: size (nx + kx)
            double precision                       :: ty(7)                            ! As above, but for the y-coordinate
            double precision                       :: tz(7)                            ! As above, but for the z-coordinate
            double precision                       :: bcoef(4, 4, 4)                   ! Matrix of coefficients of the b-spline interpolant
            double precision                       :: w0(12)                           ! Work array for db2val
            double precision                       :: w1(4)                            ! Work array for db2val
            double precision                       :: w2(3, 3)                         ! Work array db3val

            double precision                       :: baseline_time1                   ! First time in the dataset for that orbit (time to reach pi/8)
            double precision                       :: baseline_time2                   ! First time in the dataset for that orbit (time to reach pi/8)
            double precision                       :: baseline_time3                   ! First time in the dataset for that orbit (time to reach pi/8)
            double precision                       :: baseline_time4                   ! First time in the dataset for that orbit (time to reach pi/8)
            double precision                       :: actual_time1                     ! Actual time in the dataset for that orbit
            double precision                       :: actual_time2                     ! Actual time in the dataset for that orbit
            double precision                       :: actual_time3                     ! Actual time in the dataset for that orbit
            double precision                       :: actual_time4                     ! Actual time in the dataset for that orbit
            double precision                       :: orbit_number(4)
            double precision                       :: nominal_time_step1               ! nominal time step for timestep <-> absolute conversions
            double precision                       :: nominal_time_step2               ! nominal time step for timestep <-> absolute conversions
            double precision                       :: nominal_time_step3               ! nominal time step for timestep <-> absolute conversions
            double precision                       :: nominal_time_step4               ! nominal time step for timestep <-> absolute conversions
            double precision                       :: t_end_idx1                       ! Index relating to the time-step in the database
            double precision                       :: t_end_idx2                       ! Index relating to the time-step in the database
            double precision                       :: t_end_idx3                       ! Index relating to the time-step in the database
            double precision                       :: t_end_idx4                       ! Index relating to the time-step in the database
            double precision                       :: t_end_abs


            logical                             :: extrap = .false.                 ! If extrapolation is OK in db2val

            ! First, initialise problem variables

            nx = 4
            ny = 4
            nz = 4
            kx = 3
            ky = 3
            kz = 3
            iknot = 0  ! Decide for me
            iflag = -1 ! Nothing has been done yet
            idx = 0
            idy = 0
            idz = 0

            ! Collect the required data
            ! Bias = 0 if we are dealing with 'interior' points

            t_end_floor = floor(t_end)
            t_end_ceil  = ceiling(t_end)

            n_mnfd_floor = floor(n_mnfd)
            n_mnfd_ceil  = ceiling(n_mnfd)

            J_floor     = floor(J) 
            J_ceil      = ceiling(J)

            ! Check if floor == ceil

            if (t_end_floor .eq. t_end_ceil) then

                t_end_ceil = t_end_floor + 1

            end if

            if (n_mnfd_floor .eq. n_mnfd_ceil) then

                n_mnfd_ceil = 1 + n_mnfd_floor

            end if

            if (J_floor .eq. J_ceil) then

                J_ceil = J_floor + 1

            end if

            bias_tend = 0 ! Since we do it by step now anyway!

            if (n_mnfd_floor .lt. 2) then ! We need to bias UP

                bias_n_mnfd = 1

            else if (n_mnfd_ceil .gt. n_mnfd_disc - 1) then ! We need to bias DOWN

                bias_n_mnfd = -1

            else

                bias_n_mnfd = 0

            end if

            if (J_floor .lt. 2) then

                bias_J = 1

            else if (J_ceil .gt. num_orbits - 1) then

                bias_J = -1

            else

                bias_J = 0

            end if

            ! Get time-step to the desired time for each orbit

            nominal_time_step1 = abs(dataset(1, n_mnfd_floor+bias_n_mnfd, J_floor+bias_J-1, 1) - &
                                     dataset(2, floor(n_mnfd), J_floor+bias_J-1, 1))
            nominal_time_step2 = abs(dataset(1, n_mnfd_floor+bias_n_mnfd, J_floor+bias_J, 1) - &
                                     dataset(2, n_mnfd_floor+bias_n_mnfd, J_floor+bias_J, 1))
            nominal_time_step3 = abs(dataset(1, n_mnfd_floor+bias_n_mnfd, J_ceil+bias_J, 1)  - &
                                     dataset(2, n_mnfd_floor+bias_n_mnfd, J_ceil+bias_J, 1))
            nominal_time_step4 = abs(dataset(1, n_mnfd_floor+bias_n_mnfd, J_ceil+bias_J+1, 1)  - &
                                     dataset(2, n_mnfd_floor+bias_n_mnfd, J_ceil+bias_J+1, 1))

            ! Get initial time to the pi/8 plane for each orbit

            baseline_time1 = dataset(1, n_mnfd_floor+bias_n_mnfd-1, J_floor+bias_J-1, 1)
            baseline_time2 = dataset(1, n_mnfd_floor+bias_n_mnfd, J_floor+bias_J, 1)
            baseline_time3 = dataset(1, n_mnfd_ceil+bias_n_mnfd, J_ceil+bias_J, 1)
            baseline_time4 = dataset(1, n_mnfd_ceil+bias_n_mnfd+1, J_ceil+bias_J+1, 1)

            ! Get the index of the array that corresponds to *approximately* the desired time for each orbit

            t_end_idx1 = floor(abs(t_end - baseline_time1) / nominal_time_step1)
            t_end_idx2 = floor(abs(t_end - baseline_time2) / nominal_time_step2)
            t_end_idx3 = ceiling(abs(t_end - baseline_time3) / nominal_time_step3)
            t_end_idx4 = ceiling(abs(t_end - baseline_time4) / nominal_time_step4)

            x_abcissae_step = (/t_end_idx1, t_end_idx2, t_end_idx3, t_end_idx4/)

            ! Use the indices above to generate what the actual time we've obtained for each orbit is (c/f floor/ceiling)

            actual_time1 = dataset(x_abcissae_step(1), n_mnfd_floor+bias_n_mnfd-1, J_floor+bias_J-1, 1)
            actual_time2 = dataset(x_abcissae_step(2), n_mnfd_floor+bias_n_mnfd, J_floor+bias_J, 1)
            actual_time3 = dataset(x_abcissae_step(3), n_mnfd_ceil+bias_n_mnfd, J_ceil+bias_J, 1)
            actual_time4 = dataset(x_abcissae_step(4), n_mnfd_ceil+bias_n_mnfd+1, J_ceil+bias_J+1, 1)

            ! Initialise abcissae

            x_abcissae = (/actual_time1, actual_time2, actual_time3, actual_time4/)
            y_abcissae      = (/ n_mnfd_floor+bias_n_mnfd - 1, n_mnfd_floor+bias_n_mnfd, n_mnfd_ceil+bias_n_mnfd, &
                                n_mnfd_ceil + bias_n_mnfd + 1/)
            z_abcissae      = (/ J_floor + bias_J - 1, J_floor+bias_J, J_ceil+bias_J, J_ceil+bias_J+1 /)

            ! The interpolation routines require abcissae to be strictly increasing. Thus, adjust the points we sample
            ! from to enforce this (note t -ve!)

            do ii = 2, size(x_abcissae)

10              if (x_abcissae(ii) .ge. x_abcissae(ii-1)) then
                    
                    x_abcissae_step(ii) = x_abcissae_step(ii) + 1
                    x_abcissae(ii) = dataset(x_abcissae_step(ii), n_mnfd_floor+bias_n_mnfd, z_abcissae(ii), 1)

                    GOTO 10

                end if

            end do

            ! The upper time-step may now not be within the range of t that we're asking for (noting that t_end is -ve!)
            ! If so, manually adjust the upper bound to be within the range

20          if (x_abcissae(4) .gt. t_end) then
                
                x_abcissae_step(4) = x_abcissae_step(4) + 1
                x_abcissae(4) = dataset(x_abcissae_step(4), n_mnfd_floor+bias_n_mnfd, z_abcissae(4), 1)

                GO TO 20

            end if

            ! Same for the lower time-step (again, t is -ve!)

30          if (x_abcissae(1) .lt. t_end) then
                
                x_abcissae_step(1) = x_abcissae_step(1) - 1
                x_abcissae(1) = dataset(x_abcissae_step(1), n_mnfd_floor+bias_n_mnfd, z_abcissae(1), 1)

                GO TO 30

            end if

            ! Make the abcissae and the time we ask for absolute to avoid strictly increasing/decreasing errors

            x_abcissae = abs(x_abcissae)
            t_end_abs = abs(t_end)

            do dimension = 2, 7 ! Loop through state vector

                ! Manually override initialisation counters for db2val since new data is computed here

                inbvx = 1
                inbvy = 1
                inbvz = 1
                iloy  = 1
                iloz  = 1

                do ii = 1, size(x_abcissae_step)

                    do jj = 1, size(y_abcissae)

                        do kk = 1, size(z_abcissae)

                            fnc_value(ii, jj, kk) = dataset(x_abcissae_step(ii), y_abcissae(jj), z_abcissae(kk), dimension)

                        end do

                    end do

                end do

                call db3ink(x_abcissae, nx, y_abcissae, ny, z_abcissae, nz, fnc_value, kx, ky, kz, &
                            iknot, tx, ty, tz, bcoef, iflag)

                if (iflag .ne. 0) then
                    
                    print *, "iflag at failure", iflag
                    error STOP "Interpolation initialisation failed."

                end if

                call db3val(t_end_abs, n_mnfd, J, idx, idy, idz, tx, ty, tz, nx, ny, nz, kx, ky, kz, bcoef, &
                            interpolated_state_vector(dimension-1), iflag, inbvx, inbvy, inbvz, iloy, &
                            iloz, w2, w1, w0, extrap) 

                if (iflag .ne. 0) then
                    
                    print *, "iflag at failure", iflag
                    error STOP "Interpolation evaluation failed."

                end if

            end do ! dimension

        end subroutine bspline_interpolate_step

        subroutine bspline_interpolate_time(t_end, n_mnfd, J, tmani)

            use variable_initialisation

            implicit none

            double precision, intent(in)           :: t_end                            ! Backwards integration time interpolant
            double precision, intent(in)           :: n_mnfd                           ! Along-orbit direction interpolant
            
            double precision, intent(in)           :: J                                ! Which orbit we're studying (bruting over J)

            double precision, intent(out)          :: tmani                            ! Final interpolated structure 

            integer                             :: bias_n_mnfd                      ! Bias to help with boundary conditions
            integer                             :: bias_tend                        ! Bias to help with boundary conditions
            integer                             :: bias_J                           ! Bias to help with boundary conditions
            integer                             :: nx, ny, nz                       ! Number of abcissae provided
            integer                             :: kx, ky, kz                       ! Order of spline pieces = 3
            integer                             :: iknot                            ! Knot sequence flag; set = 0 to automatically choose (here)
            integer                             :: iflag                            ! Status message; 0 == good
            integer                             :: t_end_floor                      ! Floor of t_end
            integer                             :: t_end_ceil                       ! Ceiling of t_end
            integer                             :: n_mnfd_floor                     ! Floor of n_mnfd
            integer                             :: n_mnfd_ceil                      ! Ceiling of n_mnfd
            integer                             :: J_floor                          ! Floor of the desired orbit
            integer                             :: J_ceil                           ! Ceiling of the desired orbit
            integer                             :: dimension                        ! Dimension of the state vector we're interpolating
            integer                             :: idx                              ! x-derivative of polynomial to evaluate: set zero to evaluate interpolant
            integer                             :: idy                              ! As above, but for y
            integer                             :: idz                              ! As above, but for z
            integer                             :: inbvx                            ! Initialisation parameter db3val
            integer                             :: inbvy                            ! Initialisation parameter db3val
            integer                             :: inbvz                            ! Initialisation parameter db3val
            integer                             :: iloy                             ! Initialisation parameter db3val
            integer                             :: iloz                             ! Initialisation parameter db3val

            double precision                       :: x_abcissae(4)                    ! Array of x-abcissae
            double precision                       :: y_abcissae(4)                    ! Array of y-abcissae
            double precision                       :: z_abcissae(4)                    ! Array of z-abcissae
            double precision                       :: fnc_value(4, 4, 4)               ! Value of function to use in interpolating
            double precision                       :: tx(7)                            ! Locations of the knots: size (nx + kx)
            double precision                       :: ty(7)                            ! As above, but for the y-coordinate
            double precision                       :: tz(7)                            ! As above, but for the z-coordinate
            double precision                       :: bcoef(4, 4, 4)                   ! Matrix of coefficients of the b-spline interpolant
            double precision                       :: w0(12)                           ! Work array for db2val
            double precision                       :: w1(4)                            ! Work array for db2val
            double precision                       :: w2(3, 3)                         ! Work array db3val

            logical                             :: extrap = .false.                 ! If extrapolation is OK in db2val

            ! First, initialise problem variables

            nx = 4
            ny = 4
            nz = 4
            kx = 3
            ky = 3
            kz = 3
            iknot = 0  ! Decide for me
            iflag = -1 ! Nothing has been done yet
            idx = 0
            idy = 0
            idz = 0

            ! Now, access the correct portions of the dataset for each dimension

            ! Collect the required data
            ! Bias = 0 if we are dealing with 'interior' points

            t_end_floor = floor(t_end)
            t_end_ceil  = ceiling(t_end)

            n_mnfd_floor = floor(n_mnfd)
            n_mnfd_ceil  = ceiling(n_mnfd)

            J_floor     = floor(J) 
            J_ceil      = ceiling(J)

            ! Check if floor == ceil

            if (t_end_floor .eq. t_end_ceil) then

                t_end_ceil = t_end_floor + 1

            end if

            if (n_mnfd_floor .eq. n_mnfd_ceil) then

                n_mnfd_ceil = 1 + n_mnfd_floor

            end if

            if (J_floor .eq. J_ceil) then

                J_ceil = J_floor + 1

            end if

            if (t_end_floor .lt. 2) then ! We need to bias UP

                bias_tend = 1

            else if (t_end_ceil .gt. t_end_disc-1) then ! We need to bias DOWN

                bias_tend = -1

            else

                bias_tend = 0

            end if

            if (n_mnfd_floor .lt. 2) then ! We need to bias UP

                bias_n_mnfd = 1

            else if (n_mnfd_ceil .gt. n_mnfd_disc - 1) then ! We need to bias DOWN

                bias_n_mnfd = -1

            else

                bias_n_mnfd = 0

            end if

            if (J_floor .lt. 2) then

                bias_J = 1

            else if (J_ceil .gt. num_orbits - 1) then

                bias_J = -1

            else

                bias_J = 0

            end if

            do dimension = 1,1 ! Loop through state vector

                ! Manually override initialisation counters for db2val since new data is computed here

                inbvx = 1
                inbvy = 1
                inbvz = 1
                iloy  = 1
                iloz  = 1

                fnc_value(:,:, :) = dataset( (t_end_floor+bias_tend-1):(t_end_ceil+bias_tend+1), &
                                          (n_mnfd_floor+bias_n_mnfd-1):(n_mnfd_ceil+bias_n_mnfd+1), &
                                            (J_floor+bias_J-1):(J_ceil+bias_J+1), dimension) ! (tend, nmnfd)

                x_abcissae(:) = (/t_end_floor+bias_tend-1, t_end_floor+bias_tend, t_end_ceil+bias_tend, t_end_ceil+bias_tend+1/)
                y_abcissae(:) = (/n_mnfd_floor+bias_n_mnfd-1, n_mnfd_floor+bias_n_mnfd, n_mnfd_ceil+bias_n_mnfd, &
                                n_mnfd_ceil+bias_n_mnfd+1/)
                z_abcissae(:) = (/J_floor+bias_J-1, J_floor+bias_J, J_ceil+bias_J, J_ceil+bias_J+1/)

                call db3ink(x_abcissae, nx, y_abcissae, ny, z_abcissae, nz, fnc_value, kx, ky, kz, iknot, tx, ty, tz, bcoef, &
                            iflag)

                if (iflag .ne. 0) then
                    
                    print *, "iflag at failure", iflag
                    error STOP "Interpolation initialisation failed."

                end if

                call db3val(t_end, n_mnfd, J, idx, idy, idz, tx, ty, tz, nx, ny, nz, kx, ky, kz, bcoef, &
                            tmani, iflag, inbvx, inbvy, inbvz, iloy, iloz, w2, w1, w0, extrap) 

                if (iflag .ne. 0) then
                    
                    print *, "iflag at failure", iflag
                    error STOP "Interpolation evaluation failed."

                end if

            end do ! dimension

        end subroutine bspline_interpolate_time

        subroutine b_spline_interpolate_perturbed_conds(n_mnfd, J, interpolated_state_vector)

            use variable_initialisation

            implicit none

            double precision, intent(in)           :: n_mnfd                           ! Along-orbit direction interpolant
            double precision, intent(in)           :: J                                ! Which orbit we're studying (bruting over J)
            double precision, intent(out)          :: interpolated_state_vector(6)     ! Final interpolated structure

            integer                             :: bias_n_mnfd                      ! Bias to help with boundary conditions
            integer                             :: bias_J                           ! Bias to help with boundary conditions
            integer                             :: nx                               ! Number of abcissae provided
            integer                             :: ny                               ! Number of abcissae provided
            integer                             :: kx                               ! Order of spline pieces = 3
            integer                             :: ky                               ! Order of spline pieces = 3
            integer                             :: iknot                            ! Knot sequence flag; set = 0 to automatically choose (here)
            integer                             :: iflag                            ! Status message; 0 == good
            integer                             :: n_mnfd_floor                     ! Floor of n_mnfd
            integer                             :: n_mnfd_ceil                      ! Ceiling of n_mnfd
            integer                             :: J_floor                          ! Floor of the desired orbit
            integer                             :: J_ceil                           ! Ceiling of the desired orbit
            integer                             :: dimension                        ! Dimension of the state vector we're interpolating
            integer                             :: idx                              ! x-derivative of polynomial to evaluate: set zero to evaluate interpolant
            integer                             :: idy                              ! As above, but for y
            integer                             :: inbvx                            ! Initialisation parameter db2val
            integer                             :: inbvy                            ! Initialisation parameter db2val
            integer                             :: iloy                             ! Initialisation parameter db2val

            double precision                       :: x_abcissae(4)                    ! Array of x-abcissae
            double precision                       :: y_abcissae(4)                    ! Array of y-abcissae
            double precision                       :: fnc_value(4, 4)                  ! Value of function to use in interpolating
            double precision                       :: tx(7)                            ! Locations of the knots: size (nx + kx)
            double precision                       :: ty(7)                            ! As above, but for the y-coordinate
            double precision                       :: bcoef(4, 4)                      ! Matrix of coefficients of the b-spline interpolant
            double precision                       :: w0(12)                           ! Work array for db2val
            double precision                       :: w1(4)                            ! Work array for db2val

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

            else if (J_ceil .gt. 99) then ! We need to bias DOWN

                bias_J = -1

            else

                bias_J = 0

            end if

            if (n_mnfd_floor .lt. 2) then ! We need to bias UP

                bias_n_mnfd = 1

            else if (n_mnfd_ceil .gt. 99) then ! We need to bias DOWN

                bias_n_mnfd = -1

            else

                bias_n_mnfd = 0

            end if

            ! Do the position vector

            do dimension = 2, 7 ! Loop through state vector

                ! Manually override initialisation counters for db2val since new data is computed here

                inbvx = 1
                inbvy = 1
                iloy  = 1

                fnc_value(:,:) = perturbed_conds_dataset((n_mnfd_floor+bias_n_mnfd-1):(n_mnfd_ceil+bias_n_mnfd+1), &
                                            (J_floor+bias_J-1):(J_ceil+bias_J+1), dimension-1) ! (tend, nmnfd)

                x_abcissae(:) = (/n_mnfd_floor+bias_n_mnfd-1, n_mnfd_floor+bias_n_mnfd, n_mnfd_ceil+bias_n_mnfd, &
                                n_mnfd_ceil+bias_n_mnfd+1/)
                y_abcissae(:) = (/J_floor+bias_J-1, J_floor+bias_J, J_ceil+bias_J, J_ceil+bias_J+1/)

                call db2ink(x_abcissae, nx, y_abcissae, ny, fnc_value, kx, ky, iknot, tx, ty, bcoef, iflag)

                if (iflag .ne. 0) error STOP "Interpolation initialisation failed."

                call db2val(n_mnfd, J, idx, idy, tx, ty, nx, ny, kx, ky, bcoef, interpolated_state_vector(dimension-1), iflag, &
                            inbvx, inbvy, iloy, w1, w0, extrap) 

                if (iflag .ne. 0) error STOP "Interpolation evaluation failed."
                
            end do ! dimension

            ! Do the manifold time - this is 2D!

    end subroutine b_spline_interpolate_perturbed_conds                

end module compute_spline
