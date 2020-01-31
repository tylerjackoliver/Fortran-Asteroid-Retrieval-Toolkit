module compute_spline

    use precision_kinds     ! Access to working precision types
    use bspline_module      ! Uses the spine library
    use problem_parameters  ! Grants access to the dataset

    implicit none

    contains

        subroutine bspline_interpolate(t_end, n_mnfd, J, interpolated_state_vector)

            implicit none

            integer, intent(in)                 :: t_end                            ! Backwards integration time interpolant
            integer, intent(in)                 :: n_mnfd                           ! Along-orbit direction interpolant
            
            integer, intent(in)                 :: J                                ! Which orbit we're studying (bruting over J)

            real(kind=dp), intent(out)          :: interpolated_state_vector(6)     ! Final interpolated structure 

            interpolated_state_vector = dataset(t_end, n_mnfd, J, :)

        end subroutine bspline_interpolate

end module compute_spline
