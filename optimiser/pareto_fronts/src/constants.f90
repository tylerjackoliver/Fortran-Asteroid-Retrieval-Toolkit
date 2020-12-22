! @brief Derives useful constants for use in program runtimes
module constants
  
  implicit none

  double precision, parameter    :: pi = 4.0*datan(1.0d0)                                                         ! ...pi

  double precision, parameter    :: mu = 1.32712440018e11                                                         ! GM_{Sun}; km^3/s^-2 

  double precision, parameter    :: au = 149597870.7000                                                           ! ...AU, km

  double precision, parameter    :: time_to_angle_quotient = (2.d0 * pi) / (86400.d0 * 365.26d0)                  ! Convert from a time to an angle in the synodic frame

  double precision, parameter    :: position_non_dimensionalise_quotient = 1.d0 / au                              ! Non-dimensionalise a position from km/s

  double precision, parameter    :: position_dimensionalise_quotient = au                                         ! Dimensionalise a position into km/s

  double precision, parameter    :: velocity_dimensionalise_quotient = au * 2.d0 * pi / (86400.d0 * 365.25d0)     ! Dimensionalise a velocity into km/s

  double precision, parameter    :: velocity_non_dimensionalise_quotient = 365.25d0 * 86400.d0 / (2.d0 * pi * au) ! Non-dimensionalise a velocity from km/s

  double precision, parameter    :: mu3bp = 3.0032080443d-06                                                      ! CR3BP problem parameter, Sun-Earth excl. Moon

end module constants
