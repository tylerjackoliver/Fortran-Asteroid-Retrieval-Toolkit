! ///////////////////////////////////////////////////////////////
!
! Defines parameter kinds such that a double is actually twice the 
! precision of a single, quad is actually twice the precision of a
! double.
!
! c.f. Metcalf & Reid, Modern Fortran Explained
!
! Date
! ~~~~
! Apr 2018
!
! Author
! ~~~~~~
! Self
!
! //////////////////////////////////////////////////////////////

module precision_kinds

    integer, parameter :: &

        sp = kind(1.0),                              &
        dp = selected_real_kind(2*precision(1.0_sp)),&
        qp = selected_real_kind(4*precision(1.0_sp))

end module precision_kinds
