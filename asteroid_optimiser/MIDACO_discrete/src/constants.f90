! ///////////////////////////////////////////////////////////////
!
! Defines useful constants for later calculations. Note that mu
! is NOT the CR3BP parameter, but the graviational potential of the
! Sun
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


module constants

  use precision_kinds
  
  implicit none

  real(kind=dp) :: pi = 4.0_dp*datan(1.0_dp)

  real(kind=dp) :: mu = 1.32712440018e11 ! GM_{Sun}; km^3/s^-2 

  real(kind=dp) :: au = 149597870.7000 ! km

end module constants
