module integrator

    implicit none

    interface 

        INTEGER (C_INT) function integrate_c(x, t) bind(C, name="integratorCpp")

            use iso_c_binding

            implicit none

            TYPE(C_PTR), VALUE                          :: x            ! Input/output array

            REAL(kind=C_DOUBLE)                         :: t            ! Integration time

        end function integrate_c

    end interface

    contains

        subroutine integrate(xIn, t)

            use iso_c_binding

            implicit none

            double precision, intent(inout)             :: xIn(6)
            double precision, intent(in)                :: t

            real(kind=c_double), allocatable, target    :: x(:)     ! ISO-compatible version of xIn
            real(kind=c_double)                         :: t_c

            integer                                     :: c_dum

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