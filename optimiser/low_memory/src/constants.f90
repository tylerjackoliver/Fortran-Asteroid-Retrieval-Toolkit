! ///////////////////////////////////////////////////////////////
!
! Defines useful constants for later calculations. 
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

  double precision, parameter    :: pi = 4.d0*datan(1.0d0)                                                        ! Pi

  double precision, parameter    :: mu = 1.32712440018e11                                                         ! GM_{Sun}; km^3/s^-2 

  double precision, parameter    :: au = 149597870.7000                                                           ! km

  double precision, parameter    :: time_to_angle_quotient = (2.d0 * pi) / (86400.d0 * 365.26d0)                  ! Converts from ephemeris seconds to angle swept at

  double precision, parameter    :: position_non_dimensionalise_quotient = 1.d0 / au                              ! Convert from dimensionalised position -> non-dimensional position

  double precision, parameter    :: position_dimensionalise_quotient = au                                         ! Convert from non-dimensionalised position -> dimensionalised position

  double precision, parameter    :: velocity_dimensionalise_quotient = au * 2.d0 * pi / (86400.d0 * 365.25d0)     ! Convert fromm non-dimensionalised velocity -> dimensionalised velocity

  double precision, parameter    :: velocity_non_dimensionalise_quotient = 365.25d0 * 86400.d0 / (2.d0 * pi * au) ! Convert from dimensionalised velocity to non-dimensionalised velocity

  double precision, parameter    :: mu3bp = 3.0032080443d-06                                                      ! Mass parameter in the three-body problem (Sanchez et. al., 2016b)

end module constants
