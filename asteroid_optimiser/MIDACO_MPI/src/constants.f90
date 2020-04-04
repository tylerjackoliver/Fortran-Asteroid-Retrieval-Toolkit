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

  real(kind=dp), parameter    :: pi = 4.0_dp*datan(1.0_dp)

  real(kind=dp), parameter    :: mu = 1.32712440018e11 ! GM_{Sun}; km^3/s^-2 

  real(kind=dp), parameter    :: au = 149597870.7000 ! km

  real(kind=dp), parameter    :: time_to_angle_quotient = (2.d0 * pi) / (86400.d0 * 365.26d0)

  real(kind=dp), parameter    :: position_non_dimensionalise_quotient = 1.d0 / au

  real(kind=dp), parameter    :: position_dimensionalise_quotient = au

  real(kind=dp), parameter    :: velocity_dimensionalise_quotient = au * 2.d0 * pi / (86400.d0 * 365.25d0)

  real(kind=dp), parameter    :: velocity_non_dimensionalise_quotient = 365.25d0 * 86400.d0 / (2.d0 * pi * au)

  real(kind=dp), parameter    :: mu3bp = 3.0032080443d-06

end module constants
