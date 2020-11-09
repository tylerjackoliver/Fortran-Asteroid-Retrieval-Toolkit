! -----------------------------------------------------------------
!
! INTEGRATOR MODULE
! 
! @brief Contains all the subroutines necessary to interface with a given BOOST C++ integrator
!
! ------------------------------------------------------------------
module integrator

    implicit none

    interface 

        ! @brief Calls a the C++ integrator from a C wrapping function.
        ! @param[inout] x A C-style pointer to a data array representing the state to integrate. On entry, it contains the initial state. On exit, it contains the final state.
        ! @param[in] t Integration time; assumed to begin at zero and end at t (which is negative.)
        INTEGER (C_INT) function integrate_c(x, t) bind(C, name="integratorCpp")

            use iso_c_binding

            implicit none

            TYPE(C_PTR), VALUE  :: x            ! Input/output array
            REAL(kind=C_DOUBLE) :: t            ! Integration time

        end function integrate_c

    end interface

    contains

        ! @brief Integrates a given state xIn backwards by a time t
        ! @param[inout] xIn The initial condition on entry, and the final condition on exit
        ! @param[in] t The time for which to integrate
        subroutine integrate(xIn, t)
            use iso_c_binding

            implicit none

            double precision, intent(inout)             :: xIn(6)
            double precision, intent(in)                :: t

            real(kind=c_double), allocatable, target    :: x(:)     ! ISO-compatible version of xIn
            real(kind=c_double)                         :: t_c

            integer                                     :: c_dum    ! Dummy C return variable

            allocate(x(size(xIn)))

            ! Assign
            x = xIn
            t_c = t

            ! Integrate
            c_dum = integrate_c(c_loc(x), t_c)

            ! Re-assign
            xIn = x

            ! Clean
            deallocate(x)
        end subroutine integrate

end module integrator