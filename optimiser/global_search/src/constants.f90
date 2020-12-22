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

   
  
  implicit none

  double precision, parameter    :: pi = 4.0*datan(1.0d0)

  double precision, parameter    :: mu = 1.32712440018e11 ! GM_{Sun}; km^3/s^-2 

  double precision, parameter    :: au = 149597870.7000 ! km

  double precision, parameter    :: time_to_angle_quotient = (2.d0 * pi) / (86400.d0 * 365.26d0)

  double precision, parameter    :: position_non_dimensionalise_quotient = 1.d0 / au

  double precision, parameter    :: position_dimensionalise_quotient = au

  double precision, parameter    :: velocity_dimensionalise_quotient = au * 2.d0 * pi / (86400.d0 * 365.25d0)

  double precision, parameter    :: velocity_non_dimensionalise_quotient = 365.25d0 * 86400.d0 / (2.d0 * pi * au)

  double precision, parameter    :: mu3bp = 3.0032080443d-06

end module constants
