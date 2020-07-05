!------------------------------------------------------------------------------
! Fortran Asteroid Retrieval Tool (FART) v1.0: integrator module
!------------------------------------------------------------------------------
!
! MODULE: Integrator
!
!> @author
!> Jack Tyler, University of Southampton
!
! DESCRIPTION: 
!> Contains subroutines and an interface to an integrator in ODEINT/C++
!
! REVISION HISTORY:
! 01 Apr 2020 - Initial Version
! 01 Jul 2020 - Refactoring; add Doxygen support
!------------------------------------------------------------------------------


module integrator

    implicit none

    interface 

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of interface. 
        !> @brief
        !> Defines an interface to the integratorCpp routine in integrator.cc via
        !! ISO_C_BINDING
        !
        ! REVISION HISTORY:
        ! 01 Apr 2020 - Initial version
        !
        !> @param[in] x, t
        !> @param[out] Integer status variable
        !--------------------------------------------------------------------------- 

        INTEGER (C_INT) function integrate_c(x, t) bind(C, name="integratorCpp")

            use iso_c_binding

            implicit none

            TYPE(C_PTR), VALUE                          :: x            ! Input/output array

            REAL(kind=C_DOUBLE)                         :: t            ! Integration time

        end function integrate_c

    end interface

    contains

        !---------------------------------------------------------------------------  
        !> @author 
        !> Jack Tyler, University of Southampton
        !
        ! DESCRIPTION: 
        !> Brief description of subroutine. 
        !> @brief
        !> Integrates a state xIn backwards by t. Only handles backward integration.
        !
        ! REVISION HISTORY:
        ! 01 Apr 2020 - Initial version
        !
        !> @param[in] xIn, t
        !> @param[out] xIn
        !--------------------------------------------------------------------------- 

        subroutine integrate(xIn, t)

            use iso_c_binding

            implicit none

            double precision, intent(inout)             :: xIn(6)   ! Input/output array
            double precision, intent(in)                :: t        ! Time to integrate backwards by

            real(kind=c_double), allocatable, target    :: x(:)     ! ISO-compatible version of xIn
            real(kind=c_double)                         :: t_c      ! ISO-compatible version of t

            integer                                     :: c_dum    ! Dummy integer status variable

            !
            ! Allocate space for the ISO-style array
            !

            allocate(x(size(xIn)))

            !
            ! Assign values over
            !

            x = xIn
            t_c = t

            !
            ! Integrate
            !

            c_dum = integrate_c(c_loc(x), t_c)      ! Pass the location of x in memory in a style C can understand

            !
            ! Re-assign
            !

            xIn = x

            !
            ! Clean
            !

            deallocate(x)

        end subroutine integrate

end module integrator